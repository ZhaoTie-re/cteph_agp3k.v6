nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// Does the association survive removing the depth difference?
//
//   down_sampling.nf levelled every platform onto the cohort's baseline depth,
//   re-called the variants, and produced one genotype fileset per replicate. This
//   pipeline runs the association on them and asks whether the result moves.
//
//   WHY IT MATTERS HERE  Platform and phenotype are the same axis in this cohort:
//   every case sits on a non-HiSeqX platform, every control on HiSeqX. So a
//   depth-driven genotyping difference would land exactly where an association
//   would, and "is this biology or is it depth?" is not a rhetorical question.
//   Take the depth difference away; if the effect is still there, it was not made
//   of depth.
//
//   FOUR ARMS  baseline + rep1/rep2/rep3. The baseline is the cohort's own
//   genotypes run through this same code — not the published number — so the arms
//   differ by genotypes alone. That it then reproduces the published result is the
//   check that the wiring is right. Three replicates because the subsampling is
//   random: one draw cannot separate a robust signal from a lucky one.
//
//   SNP-BASED   mirrors analysis/assoc_saige. Reads p.value.NA, the p WITHOUT the
//               saddlepoint correction: SPA's own convergence varies run to run,
//               and the genotypes are what is under test.
//   GENE-BASED  mirrors analysis/assoc_rvtest, MODERATE+HIGH impact stratum.
//
//   Stages:
//     00  SNP_ASSOC        SAIGE step2 per arm (null model reused, see the script)
//     01  GENE_ASSOC_PREP  per arm: the reference's own rvtest VCF with this
//                          replicate's genotypes swapped in
//     02  GENE_ASSOC_RUN   rvtest per arm x method
//     03  GT_TABLES        the genotypes each test reads, per arm (post-GT-QC)
//         ASSOC_REPORT     one log + TWO figures (SNP and gene kept apart: 970
//                          variants vs 3, and one frame would lend the thinner
//                          evidence the confidence of the thicker) + confusion
//                          matrices carrying the cohort's OWN missingness
//     04  MAKE_SUBSETS     the platform sample lists (only_* independent,
//                          minus_* nested — the report keeps them apart)
//     05  SAIGE_SUBSET     the MAIN ANALYSIS on each subset: step 1 refit + step 2.
//                          Not a crude test: random_model keeps related samples on
//                          purpose (the GRM absorbs them) and the PCs differ between
//                          cases and controls, so an unadjusted OR answers a
//                          different question and cannot be compared with the
//                          study's. ~26 min per subset, run in parallel.
//     06  PLATFORM_REPORT  forest plots + Cochran's Q / I-squared on the
//                          independent arms only
// ─────────────────────────────────────────────────────────────────────────────

params.down_sampling_dir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/down_sampling'
params.analysis_dir      = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/analysis'
params.wgs_results       = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results'

params.results_dir       = "${params.down_sampling_dir}/assoc_stability/results"
params.script_dir        = "${params.down_sampling_dir}/assoc_stability/scripts"
params.assoc_geno_dir    = "${params.down_sampling_dir}/results/08_assoc_genotypes"
params.reps              = ['rep1', 'rep2', 'rep3']
params.conda_geno        = 'cteph_geno_pro'
params.conda_saige       = 'saige'
params.saige_bin         = '/home/b/b37974/anaconda3/envs/saige/bin'
params.rvtest_bin        = '/home/b/b37974/rvtests/executable'

// ── SNP-based (SAIGE) ────────────────────────────────────────────────────────
// The null model of the reference run, reused for every arm. It is fit on
// genome-wide LD-pruned markers; the swap touches 970 of 5.15M, so refitting would
// move the null for reasons unrelated to the test.
params.saige_null       = "${params.analysis_dir}/assoc_saige/results/01.fullGrm/02.nullModel/cteph_agp3k.v6.saige.sex.10pc.fullGrm.null"
params.snp_base_plink   = "${params.wgs_results}/15_random_model_prep/refined_core/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.random_model"
params.snp_region       = 'chr16:53703963-54121941'
params.snp_chrom        = '16'
params.lead_variant     = 'chr16:53887925:T:C'

// ── Gene-based (rvtest) ──────────────────────────────────────────────────────
// The reference pipeline's own impact-filtered VCF and its own prepared
// phenotype/covariate/refFlat files. rvtest wants formats the upstream originals
// are not in, and it merely WARNS when the phenotype does not match — dropping
// every sample and writing an empty result that reads like a null finding.
params.rvtest_prep      = "${params.analysis_dir}/assoc_rvtest/results/01.rvtest_prepare"
params.rvtest_mh_vcf    = "${params.analysis_dir}/assoc_rvtest/results/03.info_filter/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.maf_lt_threshold.moderate_high.vcf.gz"
params.gene_region      = '4:76306733-76311130'   // the VCF spells chromosomes without 'chr'
params.lead_gene        = 'STBD1'
// Gene-based significance: Bonferroni over the genes the reference pipeline
// actually tested after its NumVar>=3 filter (8,566 here; 0.05/8,566 = 5.84e-6).
params.gene_sig         = 5.8e-6
params.gene_sig_label   = 'Bonferroni, 0.05/8,621 genes'
params.pheno_name       = 'pheno1'
params.covar_name       = 'sex,pc1_avg,pc2_avg,pc3_avg,pc4_avg,pc5_avg,pc6_avg,pc7_avg,pc8_avg,pc9_avg,pc10_avg'
// ── Platform stratification ──────────────────────────────────────────────────
// Levelling the depth does not level the platform: a T7 read at baseline depth is
// still a 100 bp DNBSEQ read. Since every control is HiSeqX and no case is, any
// non-depth platform difference is confounded with the phenotype exactly as depth
// was — so the depth result on its own does not clear the signal.
params.cram_info        = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tmp/cram.v6/cram.v6.summary.csv'
params.platform_meta    = "${params.down_sampling_dir}/tuning.fraction/results/01_platform_fractions/platform_fractions.tsv"
params.sample_id_col    = 'ID_JHRPv6'
params.platform_col     = 'WGS_Platform'
// step 1 inputs, taken from the reference run so the model is built the same way
// rather than a lookalike of it.
params.saige_prune_in   = "${params.analysis_dir}/assoc_saige/results/00.prep/01.ldPruning/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.random_model.ldPruning.prune.in"
params.saige_merge_py   = "${params.analysis_dir}/assoc_saige/scripts/merge_pheno_cov.py"
params.saige_pheno_file = "${params.wgs_results}/16_cov_pheno_prep/refined_core/popgmm_subset_on_bbj_pcs.pheno.tsv"
params.saige_cov_file   = "${params.wgs_results}/16_cov_pheno_prep/refined_core/popgmm_subset_on_bbj_pcs.cov.sex.tsv"
params.saige_pheno_col  = 'PHENO1'
params.saige_covar      = 'SEX,PC1_AVG,PC2_AVG,PC3_AVG,PC4_AVG,PC5_AVG,PC6_AVG,PC7_AVG,PC8_AVG,PC9_AVG,PC10_AVG'
params.saige_step1_threads = 32

params.gene_methods     = [
    [ tag: 'skato',   opt: '--kernel skato'   ],
    [ tag: 'cmc',     opt: '--burden cmc'     ],
    [ tag: 'zeggini', opt: '--burden zeggini' ],
]

// ─────────────────────────────────────────────────────────────────────────────
process SNP_ASSOC {

    tag "${arm}"

    executor 'slurm'
    queue    'gr10478b'
    time     '12h'

    publishDir "${params.results_dir}/00_snp_assoc", mode: 'symlink'

    input:
    tuple val(arm), path(bed), path(bim), path(fam)
    path script

    output:
    tuple val(arm), path("${arm}.snp.assoc.txt"), emit: assoc

    script:
    """
    export PATH=${params.saige_bin}:/home/b/b37974/:\$PATH
    source activate ${params.conda_saige}
    set -euo pipefail

    bash ${script} \\
        ${arm} \\
        ${bed.baseName} \\
        ${params.saige_null} \\
        ${params.snp_chrom} \\
        ${arm}.snp.assoc.txt
    """
}

// ─────────────────────────────────────────────────────────────────────────────
process GENE_ASSOC_PREP {

    tag "${arm}"

    executor 'slurm'
    queue    'gr10478b'
    time     '4h'

    publishDir "${params.results_dir}/01_gene_assoc_vcf", mode: 'symlink'

    input:
    tuple val(arm), val(geno_prefix), path(geno_files), path(ds_list)
    path script

    output:
    tuple val(arm), path("${arm}.gene.vcf.gz"), path("${arm}.gene.vcf.gz.tbi"), emit: vcf
    path "${arm}.ac.tsv", emit: ac

    script:
    def geno = geno_prefix ?: 'NONE'
    def ds   = ds_list.name == 'NO_FILE' ? 'NONE' : ds_list.name
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_geno}
    set -euo pipefail

    bash ${script} \\
        ${arm} \\
        ${params.rvtest_mh_vcf} \\
        ${params.gene_region} \\
        ${geno} \\
        ${ds} \\
        ${arm}.gene.vcf.gz
    """
}

// ─────────────────────────────────────────────────────────────────────────────
process GENE_ASSOC_RUN {

    tag "${arm} | ${method.tag}"

    executor 'slurm'
    queue    'gr10478b'
    time     '12h'

    publishDir "${params.results_dir}/02_gene_assoc", mode: 'symlink'

    input:
    tuple val(arm), path(vcf), path(tbi), val(method)
    path script

    output:
    tuple val(arm), val(method.tag), path("${arm}.${method.tag}.*.assoc"), emit: assoc
    path "${arm}.${method.tag}.log"

    script:
    """
    export PATH=${params.rvtest_bin}:/home/b/b37974/:\$PATH
    source activate ${params.conda_geno}
    set -euo pipefail

    bash ${script} \\
        ${arm} \\
        ${method.tag} \\
        '${method.opt}' \\
        ${vcf} \\
        ${params.rvtest_prep}/popgmm_subset_on_bbj_pcs.pheno.pheno_rvt.tsv \\
        ${params.pheno_name} \\
        ${params.rvtest_prep}/popgmm_subset_on_bbj_pcs.cov.sex.covar_rvt.tsv \\
        ${params.covar_name} \\
        ${params.rvtest_prep}/refFlat.hg38.nochr.txt.gz \\
        ${arm}.${method.tag}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: the genotypes the two tests actually read, per arm, so the report can
//   show what the down-sampling did to them and not only what it did to the
//   p-value. Both sides post-GT-QC — that is what the associations ran on.
// ─────────────────────────────────────────────────────────────────────────────
process GT_TABLES {

    tag "${arm}"

    executor 'slurm'
    queue    'gr10478b'
    time     '2h'

    publishDir "${params.results_dir}/03_genotypes", mode: 'symlink'

    input:
    tuple val(arm), path(bed), path(bim), path(fam), path(gene_vcf), path(gene_tbi)
    path script

    output:
    tuple val(arm), path("${arm}.snp_gt.tsv"), path("${arm}.gene_gt.tsv"), emit: gt

    script:
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_geno}
    set -euo pipefail

    bash ${script} \\
        ${arm} \\
        ${bed.baseName} \\
        ${params.lead_variant} \\
        ${gene_vcf} \\
        ${arm}.snp_gt.tsv \\
        ${arm}.gene_gt.tsv
    """
}

// ─────────────────────────────────────────────────────────────────────────────
process ASSOC_REPORT {

    executor 'slurm'
    queue    'gr10478b'
    time     '2h'

    publishDir "${params.results_dir}/03_report", mode: 'symlink'

    input:
    path snp_files
    path gene_files
    path ac_files
    path gt_files
    val  snp_spec
    val  gene_spec
    val  snp_gt_spec
    val  gene_gt_spec
    path script

    output:
    path "assoc_stability.log",           emit: log
    path "assoc_stability.snp.png",       emit: fig_snp
    path "assoc_stability.gene.png",      emit: fig_gene
    path "assoc_stability.snp.tsv",       emit: tsv_snp
    path "assoc_stability.gene.tsv",      emit: tsv_gene
    path "assoc_stability.confusion.tsv", emit: tsv_conf

    script:
    def snp_args  = snp_spec.collect     { s -> "--snp ${s}" }.join(' ')
    def gene_args = gene_spec.collect    { s -> "--gene ${s}" }.join(' ')
    def ac_args   = ac_files.collect     { f -> "--ac ${f}" }.join(' ')
    def sgt_args  = snp_gt_spec.collect  { s -> "--snp-gt ${s}" }.join(' ')
    def ggt_args  = gene_gt_spec.collect { s -> "--gene-gt ${s}" }.join(' ')
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_geno}
    set -euo pipefail

    python3 ${script} \\
        ${snp_args} \\
        ${gene_args} \\
        ${ac_args} \\
        ${sgt_args} \\
        ${ggt_args} \\
        --lead-variant ${params.lead_variant} \\
        --lead-gene ${params.lead_gene} \\
        --gene-sig ${params.gene_sig} \\
        --gene-sig-label "${params.gene_sig_label}" \\
        --out-log assoc_stability.log \\
        --out-fig-snp assoc_stability.snp.png \\
        --out-fig-gene assoc_stability.gene.png \\
        --out-tsv-snp assoc_stability.snp.tsv \\
        --out-tsv-gene assoc_stability.gene.tsv \\
        --out-tsv-confusion assoc_stability.confusion.tsv
    """
}

// ─────────────────────────────────────────────────────────────────────────────
process BUILD_SNP_BASELINE {

    tag "baseline"

    executor 'slurm'
    queue    'gr10478b'
    time     '2h'

    publishDir "${params.results_dir}/00_snp_baseline", mode: 'symlink'

    output:
    tuple val('baseline'), path("baseline.bed"), path("baseline.bim"), path("baseline.fam"), emit: plink

    script:
    def r = params.snp_region
    def chr = r.tokenize(':')[0].replace('chr', '')
    def span = r.tokenize(':')[1]
    def s = span.tokenize('-')[0]
    def e = span.tokenize('-')[1]
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_geno}
    set -euo pipefail

    # The cohort's own genotypes, cut to exactly the variants the replicates carry,
    # so the arms differ by genotypes and nothing else.
    awk -F'\\t' '\$1==${chr} && \$4>=${s} && \$4<=${e} {print \$2}' \\
        ${params.snp_base_plink}.bim > vars.txt
    n=\$(wc -l < vars.txt)
    if [ "\$n" -eq 0 ]; then
        echo "[baseline] no variant of ${params.snp_base_plink}.bim lies in ${r}" >&2
        exit 1
    fi
    plink2 --bfile ${params.snp_base_plink} --extract vars.txt --make-bed --out baseline
    echo "[baseline] \$n variants x \$(wc -l < baseline.fam) samples" >&2
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: the platform sample lists. only_* share no case (independent); minus_*
//   overlap heavily (nested). The manifest records which is which so the report
//   cannot treat a nested arm as an independent confirmation.
// ─────────────────────────────────────────────────────────────────────────────
process MAKE_SUBSETS {

    executor 'slurm'
    queue    'gr10478b'
    time     '1h'

    publishDir "${params.results_dir}/04_subsets", mode: 'symlink'

    input:
    tuple val(arm), path(bed), path(bim), path(fam)
    path script

    output:
    path "*.keep",              emit: keeps
    path "subset_manifest.tsv", emit: manifest

    script:
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_geno}
    set -euo pipefail

    python3 ${script} \\
        --fam ${fam} \\
        --cram-info ${params.cram_info} \\
        --sample-id-col "${params.sample_id_col}" \\
        --platform-col "${params.platform_col}" \\
        --out-dir . \\
        --out-manifest subset_manifest.tsv
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: the main analysis on ONE subset — SAIGE step 1 refit + step 2.
// ─────────────────────────────────────────────────────────────────────────────
process SAIGE_SUBSET {

    tag "${keep.baseName}"

    executor 'slurm'
    queue    'gr10478b'
    time     '12h'

    publishDir "${params.results_dir}/05_subset_assoc", mode: 'symlink'

    input:
    path keep
    tuple val(arm), path(bed), path(bim), path(fam)
    path script

    output:
    tuple val("${keep.baseName}"), path("${keep.baseName}.assoc.txt"), emit: assoc

    script:
    """
    export PATH=${params.saige_bin}:/home/b/b37974/:\$PATH
    source activate ${params.conda_saige}
    set -euo pipefail

    bash ${script} \\
        ${keep.baseName} \\
        ${keep} \\
        ${params.snp_base_plink} \\
        ${params.saige_prune_in} \\
        ${bed.baseName} \\
        ${params.saige_merge_py} \\
        ${params.saige_pheno_file} \\
        ${params.saige_cov_file} \\
        ${params.saige_pheno_col} \\
        ${params.saige_covar} \\
        ${params.snp_chrom} \\
        ${params.saige_step1_threads} \\
        ${keep.baseName}.assoc.txt
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: forest plots + heterogeneity, with the two families kept apart.
// ─────────────────────────────────────────────────────────────────────────────
process PLATFORM_REPORT {

    executor 'slurm'
    queue    'gr10478b'
    time     '2h'

    publishDir "${params.results_dir}/06_platform_report", mode: 'symlink'

    input:
    path assoc_files
    path manifest
    val  specs
    path script

    output:
    path "platform_stratified.log", emit: log
    path "platform_stratified.tsv", emit: tsv
    path "platform_stratified.png", emit: figure

    script:
    def args = specs.collect { x -> "--assoc ${x}" }.join(' ')
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_geno}
    set -euo pipefail

    python3 ${script} \\
        ${args} \\
        --manifest ${manifest} \\
        --platform-meta ${params.platform_meta} \\
        --variant ${params.lead_variant} \\
        --out-log platform_stratified.log \\
        --out-tsv platform_stratified.tsv \\
        --out-fig platform_stratified.png
    """
}

// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // An undefined params.x is not an error in Nextflow — it renders as the string
    // 'null' and the task fails much later, inside plink or SAIGE, with a message
    // about the tool's arguments rather than about the missing parameter. A typo in
    // a params name cost eight 26-minute jobs to discover once; check it here, where
    // it costs two seconds.
    [ 'saige_prune_in', 'saige_merge_py', 'saige_pheno_file', 'saige_cov_file',
      'saige_pheno_col', 'saige_covar', 'saige_step1_threads', 'saige_null',
      'snp_base_plink', 'snp_region', 'snp_chrom', 'lead_variant',
      'cram_info', 'platform_meta', 'rvtest_mh_vcf', 'rvtest_prep',
      'gene_region', 'lead_gene', 'pheno_name', 'covar_name' ].each { k ->
        if (params[k] == null || "${params[k]}".trim() in ['', 'null']) {
            error "params.${k} is not set. Nextflow would render it as the literal " +
                  "string 'null' and the failure would surface inside a tool, far from here."
        }
    }
    // Same for the files those params point at: a path that does not exist fails at
    // staging time with the same kind of distant, misattributed error.
    [ 'saige_prune_in', 'saige_merge_py', 'saige_pheno_file', 'saige_cov_file',
      'platform_meta', 'cram_info', 'rvtest_mh_vcf' ].each { k ->
        file(params[k], checkIfExists: true)
    }
    file("${params.saige_null}.rda",                checkIfExists: true)
    file("${params.saige_null}.varianceRatio.txt",  checkIfExists: true)
    file("${params.snp_base_plink}.bed",            checkIfExists: true)


    snp_script   = file("${params.script_dir}/snp_assoc.sh",              checkIfExists: true)
    prep_script  = file("${params.script_dir}/gene_assoc_prep.sh",        checkIfExists: true)
    run_script   = file("${params.script_dir}/gene_assoc_run.sh",         checkIfExists: true)
    rep_script   = file("${params.script_dir}/assoc_stability_report.py", checkIfExists: true)
    gt_script    = file("${params.script_dir}/gt_tables.sh",               checkIfExists: true)
    subs_script  = file("${params.script_dir}/make_platform_subsets.py",   checkIfExists: true)
    saige_script = file("${params.script_dir}/saige_subset.sh",            checkIfExists: true)
    pfrep_script = file("${params.script_dir}/platform_report.py",         checkIfExists: true)
    // A real file, because Nextflow stages path inputs and the baseline arm has no
    // genotypes to swap: it IS the reference.
    no_file      = file("${params.script_dir}/NO_FILE", checkIfExists: true)

    // ── SNP: the cohort's own genotypes, plus one arm per replicate ──────────
    baseline_ch = BUILD_SNP_BASELINE().plink
    snp_reps_ch = channel.fromList(params.reps)
        .map { rep ->
            def p = "${params.assoc_geno_dir}/${rep}.snp_based"
            tuple(rep, file("${p}.bed", checkIfExists: true),
                       file("${p}.bim", checkIfExists: true),
                       file("${p}.fam", checkIfExists: true))
        }
    snp_ch = SNP_ASSOC(baseline_ch.mix(snp_reps_ch), snp_script).assoc

    // ── Platform: the main analysis re-run on each platform subset ──────────
    subs = MAKE_SUBSETS(baseline_ch, subs_script)
    sub_ch = SAIGE_SUBSET(subs.keeps.flatten(), baseline_ch.first(), saige_script).assoc
    PLATFORM_REPORT(sub_ch.map { _t, f -> f }.collect(),
                    subs.manifest,
                    sub_ch.map { t, f -> "${t}:${f.name}" }.collect(),
                    pfrep_script)

    // ── Gene ────────────────────────────────────────────────────────────────
    gene_reps_ch = channel.fromList(params.reps)
        .map { rep ->
            def p = "${params.assoc_geno_dir}/${rep}.gene_based"
            tuple(rep, p,
                  [file("${p}.bed", checkIfExists: true),
                   file("${p}.bim", checkIfExists: true),
                   file("${p}.fam", checkIfExists: true)],
                  file("${p}.ds_samples.txt", checkIfExists: true))
        }
    gene_base_ch = channel.of( tuple('baseline', 'NONE', [], no_file) )
    prep_out = GENE_ASSOC_PREP(gene_base_ch.mix(gene_reps_ch), prep_script)

    gene_ch = GENE_ASSOC_RUN(prep_out.vcf.combine(channel.fromList(params.gene_methods)),
                             run_script).assoc

    // ── the genotypes behind both tests, per arm ────────────────────────────
    // baseline's SNP plink comes from BUILD_SNP_BASELINE, the replicates' from
    // 08_assoc_genotypes; the gene VCFs all come from GENE_ASSOC_PREP. Join on arm.
    gt_in_ch = baseline_ch.mix(snp_reps_ch)
        .join(prep_out.vcf)
        .map { arm, bed, bim, fam, vcf, tbi -> tuple(arm, bed, bim, fam, vcf, tbi) }
    gt_ch = GT_TABLES(gt_in_ch, gt_script).gt

    // ── report ──────────────────────────────────────────────────────────────
    ASSOC_REPORT(snp_ch.map  { _a, f -> f }.collect(),
                 gene_ch.map { _a, _m, f -> f }.collect(),
                 prep_out.ac.collect(),
                 gt_ch.flatMap { _a, s1, g1 -> [s1, g1] }.collect(),
                 snp_ch.map  { arm, f -> "${arm}:${f.name}" }.collect(),
                 gene_ch.map { arm, m, f -> "${arm}:${m}:${f.name}" }.collect(),
                 gt_ch.map { arm, s1, _g -> "${arm}:${s1.name}" }.collect(),
                 gt_ch.map { arm, _s, g1 -> "${arm}:${g1.name}" }.collect(),
                 rep_script)
}
