#!/usr/bin/env nextflow
/* ============================================================================
 * assoc_hla — HLA allele and amino-acid association
 *
 * WHAT THIS IS. hla.typing calls 2-field HLA alleles and IMGT-numbered amino-acid
 * residues for every sample. This tests them against CTEPH case/control status,
 * over the same three nested cohorts the SNP scans use, under BOTH a fixed-effects
 * model and a full-GRM mixed model, at BOTH the allele and the residue level.
 *
 * THE ONE IDEA THAT MAKES IT CHEAP. The dosages are exact integer counts, so they
 * are written as biallelic markers on real GRCh38 chr6 coordinates (REF = A,
 * absent; ALT = P, present — the SNP2HLA convention). plink2 and SAIGE then read
 * them like any other genotypes, and everything downstream speaks the canonical
 * 16-column plink2 --glm schema the sibling components already use.
 *
 * Placing the markers genuinely on chr6 buys one real thing beyond plotting:
 * SAIGE's --LOCO=TRUE holds chr6 out of the null they are tested against, so the
 * GRM cannot contain the signal it is meant to control for.
 *
 * METHOD. SNP2HLA / Okada / Hirata, as codified in Nat. Protoc. 2023
 * (doi:10.1038/s41596-023-00853-4), which names PLINK and SAIGE as the engines.
 * Four layers: single-allele, single-residue, the amino-acid OMNIBUS at each
 * position, and forward stepwise conditioning. See docs/METHODS.md.
 *
 * PHENOTYPE CODING. pheno_covar.tsv carries 0/1, which is what SAIGE wants.
 * plink2's own convention is 1 = control, 2 = case, and it reads 0 as MISSING --
 * so without --1 it silently finds 419 controls and no cases at all. --1 tells it
 * the column is 0/1. The two engines disagreeing about this is exactly the kind of
 * thing that produces a clean-looking run of nonsense.
 *
 * Resources go through `--rsc` only. Tools read `ext.threads`: `task.cpus` is
 * unusable here because it defaults to 1, so an elvis fallback silently runs
 * single-threaded.
 * ==========================================================================*/

nextflow.enable.dsl = 2

// =============================================================================
// SITE CONFIGURATION
// =============================================================================

// ---- Paths ----
params.project_dir       = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6'
params.base_dir          = "${params.project_dir}/analysis/assoc_hla"
params.script_dir        = "${params.base_dir}/scripts"
params.shared_script_dir = "${params.project_dir}/analysis/_shared/scripts"
params.out_dir           = "${params.base_dir}/results"
params.model_inputs      = "${params.project_dir}/wgs.auto.par/results/12_model_inputs"
params.hla_dir           = "${params.project_dir}/hla.typing/results"
params.info_dir          = "${params.project_dir}/info"

params.plink2        = '/home/b/b37974/plink2'
params.tool_dir      = file(params.plink2).parent.toString()
params.conda_env     = 'cteph_geno_pro'
params.rscript       = '/home/b/b37974/anaconda3/envs/r_work/bin/Rscript'
params.saige_rscript = '/home/b/b37974/anaconda3/envs/saige/bin/Rscript'
params.saige_step1   = '/home/b/b37974/anaconda3/envs/saige/bin/step1_fitNULLGLMM.R'
params.saige_step2   = '/home/b/b37974/anaconda3/envs/saige/bin/step2_SPAtests.R'

// ---- Design ----
// The same three nested cohorts as every other association component, narrowest
// first. Every one of their samples is HLA-typed — 2,195 / 2,508 / 3,101, i.e.
// 100 % — so nothing is lost to the join and these results sit on exactly the
// sample sets the genome-wide scans used.
params.Cohorts = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']

// The two estimators. This is the comparison axis: same samples, same markers,
// same covariates, the GRM present or absent.
//   fixed  — plink2 --glm, relatives REMOVED (the fixed_model sample set), the
//            project's existing fixed-effects convention
//   random — SAIGE step 1/2, relatives KEPT and absorbed by a full GRM
params.Models = ['fixed', 'random']
params.MarkerClasses = ['allele', 'residue']

params.NPcs        = 10
params.PcLabel     = 'bbj_mainland'
params.CovarFile   = "${params.PcLabel}_pc.sex.tsv"
params.CovarNames  = (['SEX'] + (1..params.NPcs).collect { "PC${it}_AVG" }).join(',')
params.CovarLabel  = "SEX + ${params.NPcs} ${params.PcLabel} PCs"
params.ModelLabels = ['fixed' : "${params.CovarLabel}, relatives removed",
                      'random': "${params.CovarLabel}, full GRM"]
// The long form names the whole model and belongs in a caption. As an AXIS label
// it is wider than the panel and overflows the figure, so the COMPARISON figures
// take a short form; the per-cohort figures keep the long one in their titles.
params.ModelLabelsShort = ['fixed': 'fixed effects', 'random': 'full GRM']
params.PhenoName   = 'PHENO1'
params.SampleIdCol = 'IID'
params.FixedModelDir = 'fixed_model'   // the relatedness-pruned sample set

// ---- Which markers are tested ----
// DRB3/DRB4/DRB5 are held out of the primary analysis for three independent
// reasons, none of them a judgement call:
//   1. hla.typing OPEN_QUESTIONS §1 — hemizygotes are encoded as homozygotes, so
//      their dosages read 2 where they should read 1. 305 of 2,544 allele columns
//      and 635 of 3,328 residue columns.
//   2. §2 — DRB3 fails the jMorp frequency reconciliation outright (r ≈ 0.64).
//   3. DRB3 and DRB4 have NO GRCh38 primary-assembly locus. They exist only on
//      ALT haplotypes, so any coordinate given to them would be invented.
// They still run, as a separately labelled sensitivity arm; see params.Arms.
params.DropGenes = 'DRB3,DRB4,DRB5'
params.AllGenes  = 'A,B,C,DMA,DMB,DOA,DOB,DPA1,DPB1,DQA1,DQB1,DRA,DRB1,DRB3,DRB4,DRB5,E,F,G'

// The frequency floor is a COUNT. Nat. Protoc. 2023 recommends MAF >= 1 %, but
// its reason is imputation noise and these calls are typed directly. What binds
// here is the case count, 419–439: MAF 1 % is ~44 chromosomes and about 8 on the
// case side, too thin to estimate an odds ratio from.
params.MinMac        = 20
params.MinMaf        = 0.0
params.MinCallRate   = 0.95
params.MinDetermined = 0.95   // per-position residue determination, for the omnibus

// ---- Prior literature ----
// The only published HLA association for CTEPH in a Japanese population, and the
// stratum it was found in. Kominami et al. J Hum Genet 2009;54:108-114 report
// HLA-B*52:01 (OR 2.47) and HLA-DPB1*02:02 (OR 5.07) in 99 CTEPH patients WITHOUT
// deep vein thrombosis against 380 controls, and explicitly NOT in the 61 with it.
// Screening those claims is a separate question from the scan and carries its own,
// much smaller, multiple-testing burden -- two pre-specified alleles, a stated
// direction and a stated magnitude.
params.SampleSheet      = "${params.info_dir}/cteph_agp3k.v6.20260507.xlsx"
params.SampleSheetTab   = 'Sheet1'
params.SampleSheetIdCol = 'WGS_ID'   // the ONLY id column that matches the .fam
params.PriorReports     = "${params.base_dir}/info/prior_reports.tsv"

// ---- Significance ----
// Three thresholds, always reported together. The primary call is the effective
// number of independent tests (Li & Ji 2005): 5e-8 is calibrated for a million
// independent variants and is mis-calibrated — not conservative — inside 3.4 Mb
// of the most extreme LD in the genome, while plain Bonferroni pretends ~1,000
// near-collinear markers are independent. See scripts/effective_tests.py.
params.Alpha        = 0.05
params.PGenomeWide  = 5e-8
params.MaxCondRounds = 5
params.PValue       = 'spa'   // SAIGE: 'spa' | 'no_spa'

// ---- GRM ----
params.HighLdFile = "${params.info_dir}/high-LD-regions-hg38-GRCh38_modified.txt"
params.PruneWindow = 500
params.PruneStep   = 50
params.PruneR2     = 0.2
params.VarianceRatioMarkers = 200

def cohortInput(cohort) {
    def gdir = file("${params.model_inputs}/${cohort}/genotype/random_model")
    def bed  = gdir.list().findAll { n -> n.endsWith('.bed') }
    if (!bed) error "No random_model .bed under ${gdir} for cohort ${cohort}"
    def pfx = bed[0] - '.bed'
    def fdir = file("${params.model_inputs}/${cohort}/genotype/${params.FixedModelDir}")
    def ffam = fdir.exists() ? fdir.list().findAll { n -> n.endsWith('.fam') }.sort() : []
    return [cohort,
            file("${gdir}/${pfx}.bed"), file("${gdir}/${pfx}.bim"), file("${gdir}/${pfx}.fam"),
            file("${params.model_inputs}/${cohort}/phenotype/pheno.tsv"),
            file("${params.model_inputs}/${cohort}/covariates/${params.CovarFile}"),
            ffam ? file("${fdir}/${ffam[0]}") : file('NO_FILE')]
}

def script_file(name) { file("${params.script_dir}/${name}", checkIfExists: true) }
def shared_file(name) { file("${params.shared_script_dir}/${name}", checkIfExists: true) }

// =============================================================================
// PROCESSES
// =============================================================================

/* STEP 1 · GENE_COORDS -> 00.prep/ — GRCh38 coordinates for the HLA genes, from
 * Ensembl rather than typed in by hand, because every marker position derives
 * from them. DRB3 and DRB4 drop out here: they have no primary-assembly locus. */
process GENE_COORDS {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      'coords'
    publishDir "${params.out_dir}/00.prep", mode: 'copy'

    input:
    path script

    output:
    path('hla_gene_coords.tsv'), emit: coords

    script:
    """
    set -euo pipefail
    ${params.rscript} ${script} --genes ${params.AllGenes} --out hla_gene_coords.tsv
    """
}

/* STEP 2 · BUILD_MARKERS -> <cohort>/00.prep/ — the dosage matrices become a
 * chr6 PLINK fileset. The VCF is the intermediate on purpose: it states REF and
 * ALT explicitly, so the effect allele is declared rather than inferred from a
 * bit-packing convention, and plink2 writes the .bed itself. */
process BUILD_MARKERS {
    executor 'slurm'
    queue    'gr10478b'
    time     '4h'
    tag      "${cohort}"
    publishDir { "${params.out_dir}/${cohort}/00.prep" }, mode: 'copy', pattern: 'marker_map.tsv'

    input:
    tuple val(cohort), path(fam), path(allele_dosage), path(residue_dosage),
          path(pgroup_map), path(coords)
    path script
    path residue_reference

    output:
    tuple val(cohort), path('hla_markers.bed'), path('hla_markers.bim'),
          path('hla_markers.fam'), path('marker_map.tsv'), emit: markers

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    export PATH="${params.tool_dir}:\$PATH"
    python3 ${script} \\
        --allele-dosage ${allele_dosage} --residue-dosage ${residue_dosage} \\
        --gene-coords ${coords} --fam ${fam} --pgroup-map ${pgroup_map} \\
        --residue-reference ${residue_reference} \\
        --drop-genes '${params.DropGenes}' --cohort ${cohort} \\
        --out-vcf hla_markers.vcf --out-map marker_map.tsv
    ${params.plink2} --vcf hla_markers.vcf --double-id --make-bed \\
        --out hla_markers --threads ${task.ext.threads ?: 1}
    # The .fam from --vcf carries no phenotype or sex; the association reads both
    # from --pheno/--covar, so this only has to agree on the sample set.
    test \$(wc -l < hla_markers.fam) -eq \$(wc -l < ${fam})
    """
}

/* STEP 3 · MARKER_QC -> <cohort>/00.prep/ — which markers are tested, and why
 * every one that is not was dropped. One --extract list per marker class, so the
 * allele scan and the residue scan carry separate multiple-testing burdens. */
process MARKER_QC {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      "${cohort}"
    publishDir { "${params.out_dir}/${cohort}/00.prep" }, mode: 'copy', pattern: 'marker_qc.tsv'

    input:
    tuple val(cohort), path(marker_map)
    path script

    output:
    tuple val(cohort), path('marker_qc.tsv'), path('tested.allele.txt'),
          path('tested.residue.txt'), emit: qc

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} --marker-map ${marker_map} --cohort ${cohort} \\
        --min-mac ${params.MinMac} --min-maf ${params.MinMaf} \\
        --min-call-rate ${params.MinCallRate} --min-determined ${params.MinDetermined} \\
        --arm primary --out-qc marker_qc.tsv --out-prefix tested
    """
}

/* STEP 4 · PREP_PHENO_COV -> <cohort>/00.prep/ — phenotype and covariates on the
 * genotyped samples, in the genotype file's order, with the 0/1 vs 1/2 coding
 * normalised rather than assumed. */
process PREP_PHENO_COV {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      "${cohort}"
    publishDir { "${params.out_dir}/${cohort}/00.prep" }, mode: 'copy'

    input:
    tuple val(cohort), path(fam), path(pheno), path(covar)
    path script

    output:
    tuple val(cohort), path('pheno_covar.tsv'), emit: table

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} --pheno ${pheno} --covar ${covar} --fam ${fam} \\
        --pheno-col ${params.PhenoName} --id-col ${params.SampleIdCol} \\
        --out pheno_covar.tsv
    """
}

/* STEP 5 · PRUNE_MARKERS — the GRM's marker set: genome-wide, LD-pruned, high-LD
 * and inversion regions excluded, so no handful of loci dominates the
 * relatedness estimate. Not the HLA markers: the GRM describes structure, and
 * building it from the region under test is exactly the contamination
 * --LOCO exists to prevent. */
process PRUNE_MARKERS {
    executor 'slurm'
    queue    'gr10478b'
    time     '4h'
    tag      "${cohort}"
    publishDir { "${params.out_dir}/${cohort}/00.prep" }, mode: 'copy', pattern: '*.prune.in'

    input:
    tuple val(cohort), path(bed), path(bim), path(fam)

    output:
    tuple val(cohort), path('grm_markers.bed'), path('grm_markers.bim'),
          path('grm_markers.fam'), emit: markers
    path('grm.prune.in')

    script:
    """
    set -euo pipefail
    export PATH="${params.tool_dir}:\$PATH"
    awk 'BEGIN{OFS="\\t"} !/^#/ && NF>=3 {c=\$1; sub(/^chr/,"",c); print c, \$2, \$3, "highld"}' \\
        ${params.HighLdFile} > high_ld.bed
    ${params.plink2} --bed ${bed} --bim ${bim} --fam ${fam} \\
        --exclude bed1 high_ld.bed \\
        --indep-pairwise ${params.PruneWindow} ${params.PruneStep} ${params.PruneR2} \\
        --out grm --threads ${task.ext.threads ?: 1}
    ${params.plink2} --bed ${bed} --bim ${bim} --fam ${fam} --extract grm.prune.in \\
        --make-bed --out grm_markers --threads ${task.ext.threads ?: 1}
    """
}

/* STEP 6 · FIT_NULL -> <cohort>/01.assoc/random/ — SAIGE step 1, full GRM.
 * Fitted here rather than borrowed from assoc_saige so this component has no
 * cross-component runtime dependency: a sibling's refactor must not be able to
 * change what this one tests against.
 *
 * --isCateVarianceRatio=FALSE, and the reason is a measurement rather than a
 * preference. A per-MAC-category variance ratio would suit HLA markers, whose
 * frequencies span a far wider range than the markers the ratio is estimated
 * from. But SAIGE estimates it from the STEP 1 PLINK FILE, which is the
 * LD-pruned GRM marker set -- common variants, minimum MAC 51 on this cohort,
 * and ZERO markers in SAIGE's default low category of (10, 20.5]. There is
 * nothing to estimate the low category from, so asking for it would either fail
 * or silently return one category under another name. A single ratio, as
 * assoc_saige uses, is what these markers can actually support. */
process FIT_NULL {
    executor 'slurm'
    queue    'gr10478b'
    time     '24h'
    tag      "${cohort}"
    publishDir { "${params.out_dir}/${cohort}/01.assoc/random" }, mode: 'copy'

    input:
    tuple val(cohort), path(gbed), path(gbim), path(gfam), path(pheno_covar)

    output:
    tuple val(cohort), path('null.rda'), path('null.varianceRatio.txt'), emit: nullModel
    path('null.fit.log')

    script:
    """
    set -euo pipefail
    ${params.saige_rscript} ${params.saige_step1} \\
        --plinkFile=${gbed.baseName} \\
        --phenoFile=${pheno_covar} \\
        --phenoCol=${params.PhenoName} \\
        --covarColList=${params.CovarNames} \\
        --sampleIDColinphenoFile=${params.SampleIdCol} \\
        --traitType=binary \\
        --outputPrefix=null \\
        --nThreads=${task.ext.threads ?: 1} \\
        --isDiagofKinSetAsOne=True \\
        --numRandomMarkerforVarianceRatio=${params.VarianceRatioMarkers} \\
        --skipVarianceRatioEstimation=FALSE \\
        --useSparseGRMtoFitNULL=FALSE \\
        --isCateVarianceRatio=FALSE \\
        --LOCO=TRUE \\
        --IsOverwriteVarianceRatioFile=TRUE \\
        --isCovariateOffset=FALSE 2>&1 | tee null.fit.log
    """
}

/* STEP 7 · ASSOC_FIXED — plink2 --glm on the HLA markers, relatives removed.
 * firth-fallback rather than no-firth: HLA alleles are often carried by few
 * cases, and a Wald estimate at the boundary is worse than a Firth one. */
process ASSOC_FIXED {
    executor 'slurm'
    queue    'gr10478b'
    time     '2h'
    tag      "${cohort}:${mclass}"

    input:
    tuple val(cohort), val(mclass), path(bed), path(bim), path(fam), path(extract),
          path(pheno_covar), path(keep_fam)

    output:
    tuple val(cohort), val('fixed'), val(mclass), path('assoc.glm.tsv'), emit: assoc

    script:
    def keep = keep_fam.name == 'NO_FILE' ? '' : "--keep ${keep_fam}"
    """
    set -euo pipefail
    export PATH="${params.tool_dir}:\$PATH"
    ${params.plink2} --bed ${bed} --bim ${bim} --fam ${fam} \\
        --extract ${extract} ${keep} \\
        --pheno ${pheno_covar} --pheno-name ${params.PhenoName} --1 \\
        --covar ${pheno_covar} --covar-name ${params.CovarNames} \\
        --glm omit-ref hide-covar firth-fallback \\
        --ci 0.95 --threads ${task.ext.threads ?: 1} --out assoc
    cp assoc.${params.PhenoName}.glm.logistic.hybrid assoc.glm.tsv
    """
}

/* STEP 8 · ASSOC_RANDOM — SAIGE step 2 on the HLA markers. --chrom=6 with
 * --LOCO=TRUE is what holds chr6 out of the null, which is only meaningful
 * because the markers carry real chr6 coordinates. */
process ASSOC_RANDOM {
    executor 'slurm'
    queue    'gr10478b'
    time     '4h'
    tag      "${cohort}:${mclass}"

    input:
    tuple val(cohort), val(mclass), path(bed), path(bim), path(fam), path(extract),
          path(rda), path(vr)

    output:
    tuple val(cohort), val('random'), val(mclass), path('assoc.saige.txt'), emit: assoc

    script:
    """
    set -euo pipefail
    export PATH="${params.tool_dir}:\$PATH"
    # SAIGE has no --extract, so the tested subset is materialised first.
    ${params.plink2} --bed ${bed} --bim ${bim} --fam ${fam} --extract ${extract} \\
        --make-bed --out tested --threads ${task.ext.threads ?: 1}
    ${params.saige_rscript} ${params.saige_step2} \\
        --bedFile=tested.bed --bimFile=tested.bim --famFile=tested.fam \\
        --AlleleOrder=alt-first --chrom=6 \\
        --GMMATmodelFile=${rda} --varianceRatioFile=${vr} \\
        --SAIGEOutputFile=assoc.saige.txt \\
        --is_Firth_beta=TRUE --LOCO=TRUE --is_output_moreDetails=TRUE
    """
}

/* STEP 9 · TO_SUMSTATS -> <cohort>/01.assoc/<model>/<class>/ — both engines to
 * the one schema, with the effect allele asserted to be P at every marker. An
 * inverted effect allele silently inverts every odds ratio, so it is checked
 * rather than trusted. */
process TO_SUMSTATS {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      "${cohort}:${model}:${mclass}"
    publishDir { "${params.out_dir}/${cohort}/01.assoc/${model}/${mclass}" }, mode: 'copy'

    input:
    tuple val(cohort), val(model), val(mclass), path(assoc), path(bim), path(marker_map)
    path script

    output:
    tuple val(cohort), val(model), val(mclass), path('sumstats.tsv'), emit: sumstats
    path('engine_extra.tsv')

    script:
    def engine = model == 'fixed' ? 'plink2' : 'saige'
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} --assoc ${assoc} --engine ${engine} --bim ${bim} \\
        --marker-map ${marker_map} --cohort ${cohort} --model ${model} \\
        --marker-class ${mclass} --p-value ${params.PValue} \\
        --out-sumstats sumstats.tsv --out-extra engine_extra.tsv
    """
}

/* STEP 10 · EXPORT_TESTED — the tested markers as a dosage matrix, for the
 * omnibus fits and for the marker correlation the effective-test count needs. */
process EXPORT_TESTED {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      "${cohort}:${mclass}"

    input:
    tuple val(cohort), val(mclass), path(bed), path(bim), path(fam), path(extract)

    output:
    tuple val(cohort), val(mclass), path('tested.raw'), emit: raw

    script:
    """
    set -euo pipefail
    export PATH="${params.tool_dir}:\$PATH"
    ${params.plink2} --bed ${bed} --bim ${bim} --fam ${fam} --extract ${extract} \\
        --export A --out tested --threads ${task.ext.threads ?: 1}
    """
}

/* STEP 11 · OMNIBUS -> <cohort>/02.omnibus/ — the m-1 df likelihood-ratio test at
 * each amino-acid position. This is the layer that makes the analysis HLA fine
 * mapping rather than a marker scan: it asks whether the residue CONTENT of a
 * position matters, not whether one residue does. */
process OMNIBUS {
    executor 'slurm'
    queue    'gr10478b'
    time     '4h'
    tag      "${cohort}"
    publishDir { "${params.out_dir}/${cohort}/02.omnibus" }, mode: 'copy'

    input:
    tuple val(cohort), path(raw), path(pheno_covar), path(marker_qc)
    path script

    output:
    tuple val(cohort), path('omnibus.tsv'), emit: omnibus

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} --raw ${raw} --pheno-covar ${pheno_covar} \\
        --marker-qc ${marker_qc} --covars ${params.CovarNames} \\
        --pheno-name ${params.PhenoName} --sample-id-col ${params.SampleIdCol} \\
        --cohort ${cohort} --model fixed --out omnibus.tsv
    """
}

/* STEP 12 · SIGNIFICANCE -> <cohort>/03.signals/ — three thresholds, computed
 * rather than assumed, with the effective-test one marked primary. */
process SIGNIFICANCE {
    executor 'slurm'
    queue    'gr10478b'
    time     '2h'
    tag      "${cohort}:${model}:${mclass}"
    publishDir { "${params.out_dir}/${cohort}/03.signals" }, mode: 'copy',
               saveAs: { _f -> "significance.${model}.${mclass}.tsv" }

    input:
    tuple val(cohort), val(model), val(mclass), path(sumstat), path(raw)
    path script

    output:
    tuple val(cohort), val(model), val(mclass), path('significance.tsv'), emit: sig

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} --raw ${raw} --sumstat ${sumstat} --cohort ${cohort} \\
        --model ${model} --marker-class ${mclass} --alpha ${params.Alpha} \\
        --p-genomewide ${params.PGenomeWide} --out significance.tsv
    """
}

/* STEP 12b · OMNIBUS_SIGNIFICANCE -> <cohort>/03.signals/ — the position-level
 * tests need their OWN threshold. There are ~300 of them, not ~700, and they are
 * correlated differently: residues within a position are near-collinear by
 * construction, which is exactly what the omnibus collapses. Judging a position
 * against the single-marker threshold would be judging it against the wrong
 * experiment. */
process OMNIBUS_SIGNIFICANCE {
    executor 'slurm'
    queue    'gr10478b'
    time     '2h'
    tag      "${cohort}"
    publishDir { "${params.out_dir}/${cohort}/03.signals" }, mode: 'copy',
               saveAs: { _f -> 'significance.omnibus.tsv' }

    input:
    tuple val(cohort), path(omnibus), path(raw), path(marker_qc)
    path script

    output:
    tuple val(cohort), path('significance.tsv'), emit: sig

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} --mode omnibus --raw ${raw} --omnibus ${omnibus} \\
        --marker-qc ${marker_qc} --cohort ${cohort} --model fixed \\
        --marker-class residue --alpha ${params.Alpha} \\
        --p-genomewide ${params.PGenomeWide} --out significance.tsv
    """
}

/* STEP 13 · CONDITIONAL -> <cohort>/03.signals/ — forward stepwise over the whole
 * MHC. A single causal allele drags its haplotype over the threshold; this asks
 * how many SEPARATE things are going on. Conditioning is global rather than per
 * gene, because LD here spans genes. */
process CONDITIONAL {
    executor 'slurm'
    queue    'gr10478b'
    time     '8h'
    tag      "${cohort}:${model}:${mclass}"
    publishDir { "${params.out_dir}/${cohort}/03.signals" }, mode: 'copy', pattern: 'cond.*'

    input:
    tuple val(cohort), val(model), val(mclass), path(bed), path(bim), path(fam),
          path(extract), path(marker_map), path(pheno_covar), path(rda), path(vr),
          path(sig), val(threshold)
    path script
    path adapter

    output:
    tuple val(cohort), val(model), val(mclass), path('cond.*.rounds.tsv'),
          path('cond.*.signals.tsv'), emit: cond
    path('cond.*.round*.tsv'), optional: true

    script:
    def engine_args = model == 'fixed'
        ? "--plink2 ${params.plink2} --pheno ${pheno_covar} --covar ${pheno_covar} --covar-name ${params.CovarNames} --pheno-name ${params.PhenoName}"
        : "--saige-step2 ${params.saige_step2} --saige-rscript ${params.saige_rscript} --gmmat ${rda} --variance-ratio ${vr}"
    """
    set -euo pipefail
    source activate ${params.conda_env}
    export PATH="${params.tool_dir}:\$PATH"
    python3 ${script} --bfile ${bed.baseName} --extract ${extract} \\
        --marker-map ${marker_map} --cohort ${cohort} --model ${model} \\
        --marker-class ${mclass} --threshold ${threshold} \\
        --max-rounds ${params.MaxCondRounds} --adapter ${adapter} \\
        --p-value ${params.PValue} --threads ${task.ext.threads ?: 1} \\
        ${engine_args} --out-prefix cond.${model}.${mclass}
    """
}

/* STEP 14 · PLOT_HLA_SCAN -> <cohort>/figures/01.scan/ — the canonical MHC figure:
 * both marker classes on one real chr6 axis, with the gene layout below it so a
 * peak can be read off to a gene. */
process PLOT_HLA_SCAN {
    executor 'slurm'
    queue    'gr10478b'
    time     '2h'
    tag      "${cohort}:${model}"
    publishDir { "${params.out_dir}/${cohort}/figures/01.scan" }, mode: 'copy'

    input:
    // Both classes' files are called sumstats.tsv; without stageAs they collide
    // in the task directory and Nextflow refuses to stage them.
    tuple val(cohort), val(model), path(ss_allele, stageAs: 'allele.sumstats.tsv'),
          path(ss_residue, stageAs: 'residue.sumstats.tsv'), path(marker_qc),
          val(bonf_allele), val(bonf_residue), val(effective)
    path script

    output:
    path("scan.${model}.png"), emit: png
    path('*.md'), optional: true

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} --sumstat-allele ${ss_allele} --sumstat-residue ${ss_residue} \\
        --marker-qc ${marker_qc} --cohort ${cohort} --model ${model} \\
        --model-label '${params.ModelLabels[model]}' \\
        --alpha-bonferroni-allele ${bonf_allele} \\
        --alpha-bonferroni-residue ${bonf_residue} \\
        --alpha-effective ${effective} --alpha-genomewide ${params.PGenomeWide} \\
        --out-png scan.${model}.png
    """
}

/* STEP 15 · PLOT_OMNIBUS -> <cohort>/figures/02.omnibus/ — omnibus P against IMGT
 * position, per gene, with the peptide-binding domain shaded because that is
 * where a functional signal is expected to sit. */
process PLOT_OMNIBUS {
    executor 'slurm'
    queue    'gr10478b'
    time     '2h'
    tag      "${cohort}"
    publishDir { "${params.out_dir}/${cohort}/figures/02.omnibus" }, mode: 'copy'

    input:
    tuple val(cohort), path(omnibus), path(marker_qc), path(coords), val(effective)
    path script

    output:
    path('omnibus.fixed.png'), emit: png
    path('*.md'), optional: true

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    GO=\$(awk -F'\\t' 'NR>1 {printf "%s%s", sep, \$1; sep=","}' ${coords})
    python3 ${script} --omnibus ${omnibus} --marker-qc ${marker_qc} \\
        --cohort ${cohort} --model fixed \\
        --model-label '${params.ModelLabels['fixed']}' \\
        --alpha ${effective} --gene-order "\$GO" --out-png omnibus.fixed.png
    """
}

/* STEP 16 · PLOT_FREQ_QC -> <cohort>/figures/03.qc/ — the artefact check applied
 * to the hits themselves: an allele a 61,424-person Japanese panel has never seen
 * is more likely a typing error than a finding. */
process PLOT_FREQ_QC {
    executor 'slurm'
    queue    'gr10478b'
    time     '2h'
    tag      "${cohort}:${model}"
    publishDir { "${params.out_dir}/${cohort}/figures/03.qc" }, mode: 'copy'

    input:
    tuple val(cohort), val(model), path(sumstat), path(marker_qc), val(effective)
    path script

    output:
    path("freq_qc.${model}.png"), emit: png
    path('*.md'), optional: true

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} --sumstat ${sumstat} --marker-qc ${marker_qc} \\
        --cohort ${cohort} --model ${model} --alpha ${effective} \\
        --out-png freq_qc.${model}.png
    """
}

/* STEP 17 · COMPARE_MODELS -> _comparison/ — fixed against random, the axis this
 * component exists to expose: same samples, same markers, same covariates, the
 * GRM present or absent. */
process COMPARE_MODELS {
    executor 'slurm'
    queue    'gr10478b'
    time     '2h'
    tag      "${mclass}"
    publishDir "${params.out_dir}/_comparison/figures", mode: 'copy', pattern: '*.{png,md}'
    publishDir "${params.out_dir}/_comparison/tables", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(mclass), val(keys), path(sumstats, stageAs: 'ss_*.tsv')
    path script

    output:
    path("model_compare.${mclass}.png"), emit: png
    path("model_compare.${mclass}.tsv")
    path('*.md'), optional: true

    script:
    def args = [keys, sumstats].transpose().collect { k, f -> "--sumstat '" + k + "=" + f + "'" }.join(' ')
    def labels = params.Models.collect { m -> "--model-label '" + m + "=" + params.ModelLabelsShort[m] + "'" }.join(' ')
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} ${args} ${labels} --marker-class ${mclass} \\
        --cohort-order ${params.Cohorts.join(',')} \\
        --model-order ${params.Models.join(',')} \\
        --out-png model_compare.${mclass}.png --out-table model_compare.${mclass}.tsv
    """
}

/* STEP 18 · COMPARE_COHORTS -> _comparison/ — the three nested cohorts. Nested, so
 * agreement is robustness to where the ancestry boundary was drawn and NOT
 * replication; the figure says so rather than leaving it to the reader. */
process COMPARE_COHORTS {
    executor 'slurm'
    queue    'gr10478b'
    time     '2h'
    tag      "${model}:${mclass}"
    publishDir "${params.out_dir}/_comparison/figures", mode: 'copy', pattern: '*.{png,md}'
    publishDir "${params.out_dir}/_comparison/tables", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(model), val(mclass), val(cohorts), path(sumstats, stageAs: 'ss_*.tsv')
    path script

    output:
    path("cohort_compare.${model}.${mclass}.png"), emit: png
    path("cohort_compare.${model}.${mclass}.tsv")
    path('*.md'), optional: true

    script:
    def args = [cohorts, sumstats].transpose().collect { c, f -> "--sumstat '" + c + "=" + f + "'" }.join(' ')
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} ${args} --cohort-order ${params.Cohorts.join(',')} \\
        --marker-class ${mclass} --model ${model} \\
        --model-label '${params.ModelLabelsShort[model]}' \\
        --out-png cohort_compare.${model}.${mclass}.png \\
        --out-table cohort_compare.${model}.${mclass}.tsv
    """
}

/* STEP 20 · PREP_SUBGROUPS -> <cohort>/00.prep/ — DVT status for the case series,
 * from the study's own sample sheet. It is a STRATUM, not a covariate: the prior
 * Japanese CTEPH HLA report found its associations only in patients WITHOUT deep
 * vein thrombosis, so the whole-series comparison tests a different hypothesis. */
process PREP_SUBGROUPS {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      "${cohort}"
    publishDir { "${params.out_dir}/${cohort}/00.prep" }, mode: 'copy'

    input:
    tuple val(cohort), path(pheno_covar)
    path script
    path sample_sheet

    output:
    tuple val(cohort), path('subgroups.tsv'), emit: sub

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} --sample-sheet ${sample_sheet} --sheet '${params.SampleSheetTab}' \\
        --pheno-covar ${pheno_covar} --id-col ${params.SampleSheetIdCol} \\
        --sample-id-col ${params.SampleIdCol} --pheno-name ${params.PhenoName} \\
        --cohort ${cohort} --out subgroups.tsv
    """
}

/* STEP 21 · PRIOR_SCREEN -> <cohort>/04.prior/ — does the published claim hold in
 * these data? A handful of pre-specified alleles, in the stratum they were
 * reported in, with the power to detect the reported effect stated on every row
 * so that "did not replicate" can be told apart from "could not have". */
process PRIOR_SCREEN {
    executor 'slurm'
    queue    'gr10478b'
    time     '2h'
    tag      "${cohort}:${model}"
    publishDir { "${params.out_dir}/${cohort}/04.prior" }, mode: 'copy',
               saveAs: { _f -> "prior_screen.${model}.tsv" }

    input:
    tuple val(cohort), val(model), path(raw), path(pheno_covar), path(subgroups),
          path(marker_qc)
    path script
    path prior_reports

    output:
    tuple val(cohort), val(model), path('prior_screen.tsv'), emit: screen

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} --raw ${raw} --pheno-covar ${pheno_covar} \\
        --subgroups ${subgroups} --prior-reports ${prior_reports} \\
        --marker-qc ${marker_qc} --covars ${params.CovarNames} \\
        --pheno-name ${params.PhenoName} --sample-id-col ${params.SampleIdCol} \\
        --cohort ${cohort} --model ${model} --alpha 0.05 --out prior_screen.tsv
    """
}

/* STEP 19 · WRITE_RUN_MANIFEST -> _run_info/ */
process WRITE_RUN_MANIFEST {
    executor 'local'
    tag 'manifest'
    publishDir "${params.out_dir}/_run_info", mode: 'copy'

    input:
    val _trigger

    output:
    path('run_manifest.json')

    script:
    """
    cat > run_manifest.json <<'JSON'
{
  "component": "assoc_hla",
  "cohorts": ${groovy.json.JsonOutput.toJson(params.Cohorts)},
  "models": ${groovy.json.JsonOutput.toJson(params.Models)},
  "marker_classes": ${groovy.json.JsonOutput.toJson(params.MarkerClasses)},
  "covariates": "${params.CovarNames}",
  "pc_label": "${params.PcLabel}",
  "n_pcs": ${params.NPcs},
  "drop_genes": "${params.DropGenes}",
  "min_mac": ${params.MinMac},
  "min_maf": ${params.MinMaf},
  "min_call_rate": ${params.MinCallRate},
  "min_determined": ${params.MinDetermined},
  "alpha": ${params.Alpha},
  "p_genome_wide": ${params.PGenomeWide},
  "saige_p_value": "${params.PValue}",
  "max_conditional_rounds": ${params.MaxCondRounds},
  "hla_dir": "${params.hla_dir}",
  "model_inputs": "${params.model_inputs}",
  "sample_sheet": "${params.SampleSheet}",
  "prior_reports": "${params.PriorReports}"
}
JSON
    """
}

// =============================================================================
// WORKFLOW
// =============================================================================

workflow {

    // -- [1] gene coordinates ---------------------------------------------
    GENE_COORDS(script_file('hla_gene_coords.R'))

    ch_cohort = channel.fromList(params.Cohorts).map { c -> cohortInput(c) }

    // -- [2,3] markers and their QC ---------------------------------------
    ch_build = ch_cohort
        .map { c, _bed, _bim, fam, _ph, _cv, _ffam ->
               tuple(c, fam,
                     file("${params.hla_dir}/04.residues/allele_dosage.tsv", checkIfExists: true),
                     file("${params.hla_dir}/04.residues/residue_dosage.tsv", checkIfExists: true),
                     file("${params.hla_dir}/05.qc/allele_pgroup_map.tsv", checkIfExists: true)) }
        .combine(GENE_COORDS.out.coords)

    BUILD_MARKERS(ch_build, script_file('build_hla_markers.py'),
                  file("${params.hla_dir}/04.residues/residue_reference", checkIfExists: true))
    MARKER_QC(BUILD_MARKERS.out.markers.map { c, _b, _i, _f, m -> tuple(c, m) },
              script_file('hla_marker_qc.py'))

    // -- [4,5] phenotype and the GRM markers ------------------------------
    PREP_PHENO_COV(ch_cohort.map { c, _bed, _bim, fam, ph, cv, _ffam -> tuple(c, fam, ph, cv) },
                   script_file('prep_pheno_cov.py'))
    PRUNE_MARKERS(ch_cohort.map { c, bed, bim, fam, _ph, _cv, _ffam -> tuple(c, bed, bim, fam) })

    // -- [6] the null model -----------------------------------------------
    FIT_NULL(PRUNE_MARKERS.out.markers.join(PREP_PHENO_COV.out.table))

    // -- [7,8] association, one job per cohort x marker class --------------
    ch_bfile = BUILD_MARKERS.out.markers.map { c, b, i, f, _m -> tuple(c, b, i, f) }
    ch_extract = MARKER_QC.out.qc
        .flatMap { c, _qc, a, r -> [tuple(c, 'allele', a), tuple(c, 'residue', r)] }

    ch_fixed_in = ch_extract.combine(ch_bfile, by: 0)
        .map { c, mc, ex, b, i, f -> tuple(c, mc, b, i, f, ex) }
        .combine(PREP_PHENO_COV.out.table, by: 0)
        .combine(ch_cohort.map { c, _bed, _bim, _fam, _ph, _cv, ffam -> tuple(c, ffam) }, by: 0)
    ASSOC_FIXED(ch_fixed_in)

    ch_random_in = ch_extract.combine(ch_bfile, by: 0)
        .map { c, mc, ex, b, i, f -> tuple(c, mc, b, i, f, ex) }
        .combine(FIT_NULL.out.nullModel, by: 0)
    ASSOC_RANDOM(ch_random_in)

    // -- [9] one schema ----------------------------------------------------
    ch_assoc = ASSOC_FIXED.out.assoc.mix(ASSOC_RANDOM.out.assoc)
    ch_sumstat_in = ch_assoc
        .combine(BUILD_MARKERS.out.markers.map { c, _b, i, _f, m -> tuple(c, i, m) }, by: 0)
    TO_SUMSTATS(ch_sumstat_in, script_file('hla_to_sumstats.py'))

    // -- [10,11] dosages, then the omnibus ---------------------------------
    EXPORT_TESTED(ch_extract.combine(ch_bfile, by: 0)
                            .map { c, mc, ex, b, i, f -> tuple(c, mc, b, i, f, ex) })
    ch_res_raw = EXPORT_TESTED.out.raw.filter { _c, mc, _r -> mc == 'residue' }
                                      .map { c, _mc, r -> tuple(c, r) }
    OMNIBUS(ch_res_raw.combine(PREP_PHENO_COV.out.table, by: 0)
                      .combine(MARKER_QC.out.qc.map { c, qc, _a, _r -> tuple(c, qc) }, by: 0),
            script_file('hla_omnibus.py'))

    // -- [12b] the omnibus's own threshold ----------------------------------
    OMNIBUS_SIGNIFICANCE(
        OMNIBUS.out.omnibus
            .combine(ch_res_raw, by: 0)
            .combine(MARKER_QC.out.qc.map { c, qc, _a, _r -> tuple(c, qc) }, by: 0),
        script_file('effective_tests.py'))

    // -- [12] significance --------------------------------------------------
    // Keyed on cohort AND marker class: the exported dosage matrix is per class,
    // and both models of that class share it. A plain string key, because a tuple
    // key silently never emits.
    ch_sig_in = TO_SUMSTATS.out.sumstats
        .map { c, m, mc, ss -> tuple("${c}|${mc}".toString(), c, m, mc, ss) }
        .combine(EXPORT_TESTED.out.raw.map { c, mc, r -> tuple("${c}|${mc}".toString(), r) }, by: 0)
        .map { _k, c, m, mc, ss, raw -> tuple(c, m, mc, ss, raw) }
    SIGNIFICANCE(ch_sig_in, script_file('effective_tests.py'))

    // -- [13] conditional analysis -----------------------------------------
    // The threshold each run conditions against is the PRIMARY one this scan
    // computed, read out of significance.tsv rather than hard-coded, so the
    // rounds are judged by the same rule the scan was.
    ch_thr = SIGNIFICANCE.out.sig
        .map { c, m, mc, f -> tuple("${c}|${m}|${mc}".toString(), f) }
        .splitCsv(elem: 1, sep: '\t', header: true)
        .filter { _k, row -> row.primary == '1' }
        .map { k, row -> tuple(k, row.value.toString()) }

    ch_bonf = SIGNIFICANCE.out.sig
        .map { c, m, mc, f -> tuple("${c}|${m}|${mc}".toString(), f) }
        .splitCsv(elem: 1, sep: '\t', header: true)
        .filter { _k, row -> row.threshold == 'bonferroni' }
        .map { k, row -> tuple(k, row.value.toString()) }

    ch_cond_in = ch_extract
        .combine(ch_bfile, by: 0)
        .combine(BUILD_MARKERS.out.markers.map { c, _b, _i, _f, mm -> tuple(c, mm) }, by: 0)
        .combine(PREP_PHENO_COV.out.table, by: 0)
        .combine(FIT_NULL.out.nullModel, by: 0)
        .flatMap { c, mc, ex, b, i, f, mm, pc, rda, vr ->
                   params.Models.collect { m ->
                       tuple("${c}|${m}|${mc}".toString(), c, m, mc, b, i, f, ex, mm, pc, rda, vr) } }
        .combine(SIGNIFICANCE.out.sig
                     .map { c, m, mc, s -> tuple("${c}|${m}|${mc}".toString(), s) }, by: 0)
        .combine(ch_thr, by: 0)
        .map { _k, c, m, mc, b, i, f, ex, mm, pc, rda, vr, s, thr ->
               tuple(c, m, mc, b, i, f, ex, mm, pc, rda, vr, s, thr) }
    CONDITIONAL(ch_cond_in, script_file('hla_conditional.py'),
                script_file('hla_to_sumstats.py'))

    // -- [14,15,16] per-cohort figures --------------------------------------
    ch_ss = TO_SUMSTATS.out.sumstats
    ch_scan_in = ch_ss.filter { _c, _m, mc, _s -> mc == 'allele' }
        .map { c, m, _mc, s -> tuple("${c}|${m}".toString(), c, m, s) }
        .combine(ch_ss.filter { _c, _m, mc, _s -> mc == 'residue' }
                      .map { c, m, _mc, s -> tuple("${c}|${m}".toString(), s) }, by: 0)
        .map { k, c, m, sa, sr -> tuple(c, k, m, sa, sr) }
        .combine(MARKER_QC.out.qc.map { c, qc, _a, _r -> tuple(c, qc) }, by: 0)
        .map { c, _k, m, sa, sr, qc -> tuple("${c}|${m}|allele".toString(), c, m, sa, sr, qc) }
        .combine(ch_bonf, by: 0)
        .map { _k, c, m, sa, sr, qc, ba -> tuple("${c}|${m}|residue".toString(), c, m, sa, sr, qc, ba) }
        .combine(ch_bonf, by: 0)
        .map { _k, c, m, sa, sr, qc, ba, br -> tuple("${c}|${m}|allele".toString(), c, m, sa, sr, qc, ba, br) }
        .combine(ch_thr, by: 0)
        .map { _k, c, m, sa, sr, qc, ba, br, ef -> tuple(c, m, sa, sr, qc, ba, br, ef) }
    PLOT_HLA_SCAN(ch_scan_in, script_file('plot_hla_scan.py'))

    // The omnibus figure is drawn against the OMNIBUS threshold, not the
    // single-marker one: the two score different experiments.
    ch_omni_thr = OMNIBUS_SIGNIFICANCE.out.sig
        .splitCsv(elem: 1, sep: '\t', header: true)
        .filter { _c, row -> row.primary == '1' }
        .map { c, row -> tuple(c, row.value.toString()) }
    ch_omni_in = OMNIBUS.out.omnibus
        .combine(MARKER_QC.out.qc.map { c, qc, _a, _r -> tuple(c, qc) }, by: 0)
        .combine(GENE_COORDS.out.coords)
        .combine(ch_omni_thr, by: 0)
        .map { c, om, qc, co, ef -> tuple(c, om, qc, co, ef) }
    PLOT_OMNIBUS(ch_omni_in, script_file('plot_omnibus.py'))

    ch_qc_in = ch_ss.filter { _c, _m, mc, _s -> mc == 'allele' }
        .combine(MARKER_QC.out.qc.map { c, qc, _a, _r -> tuple(c, qc) }, by: 0)
        .map { c, m, _mc, s, qc -> tuple("${c}|${m}|allele".toString(), c, m, s, qc) }
        .combine(ch_thr, by: 0)
        .map { _k, c, m, s, qc, ef -> tuple(c, m, s, qc, ef) }
    PLOT_FREQ_QC(ch_qc_in, script_file('plot_hla_freq_qc.py'))

    // -- [17,18] the two comparison axes ------------------------------------
    // Barriers, deliberately: each figure needs every cohort (or every model)
    // at once, and there is nothing to pipeline.
    ch_model_cmp = ch_ss
        .map { c, m, mc, s -> tuple(mc, "${c}/${m}".toString(), s) }
        .groupTuple(by: 0)
    COMPARE_MODELS(ch_model_cmp, script_file('plot_hla_model_compare.py'))

    ch_cohort_cmp = ch_ss
        .map { c, m, mc, s -> tuple(m, mc, c, s) }
        .groupTuple(by: [0, 1])
    COMPARE_COHORTS(ch_cohort_cmp, script_file('plot_hla_cohort_compare.py'))

    // -- [20,21] the published claims, screened -----------------------------
    PREP_SUBGROUPS(PREP_PHENO_COV.out.table, script_file('prep_subgroups.py'),
                   file(params.SampleSheet, checkIfExists: true))

    ch_allele_raw = EXPORT_TESTED.out.raw
        .filter { _c, mc, _r -> mc == 'allele' }
        .map { c, _mc, r -> tuple(c, r) }
    ch_prior_in = ch_allele_raw
        .combine(PREP_PHENO_COV.out.table, by: 0)
        .combine(PREP_SUBGROUPS.out.sub, by: 0)
        .combine(MARKER_QC.out.qc.map { c, qc, _a, _r -> tuple(c, qc) }, by: 0)
        // ONE screen per cohort, not one per model. The stratified test refits on
        // a case subset and is fixed-effects by construction; a second row under
        // the 'random' label would be the same numbers wearing a name they did
        // not earn.
        .map { c, raw, pc, sg, qc -> tuple(c, 'fixed_stratified', raw, pc, sg, qc) }
    PRIOR_SCREEN(ch_prior_in, script_file('prior_screen.py'),
                 file(params.PriorReports, checkIfExists: true))

    // -- [19] provenance ----------------------------------------------------
    WRITE_RUN_MANIFEST(COMPARE_COHORTS.out.png
                          .mix(PRIOR_SCREEN.out.screen.map { _c, _m, f -> f })
                          .collect().map { _x -> 'done' })
}
