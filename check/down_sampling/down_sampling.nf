nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// Per-platform down-sampling over the target regions
//
//   Goal : the cohort was sequenced on platforms of unequal depth, and depth
//          drives what a caller can see — coverage decides which genotypes are
//          confidently called and which go missing. So any association in these
//          regions is open to the same objection: is this biology, or is it the
//          deeper platforms yielding cleaner genotypes than the shallower ones?
//
//          This pipeline answers it by removing the difference. Every platform is
//          brought down to the cohort's baseline depth, the variants are re-called
//          from those reads, and the genotypes are scored against the cohort's
//          own. Run the association on the levelled genotypes: an effect that
//          survives cannot have been manufactured by depth heterogeneity, because
//          there is none left to manufacture it.
//
//          The target regions are a parameter, not a premise. This is a test of
//          whatever signal params.regions covers — point it at another locus and
//          the same argument holds.
//
//   WHICH SAMPLES  Not "Target_Depth == 30x". The selection comes from
//          tuning.fraction's platform_fractions.tsv: every sample whose PLATFORM
//          has KEEP_FRACTION < 1 is down-sampled by ITS OWN platform fraction.
//          Platforms at KEEP_FRACTION == 1 are already at or below the baseline
//          and are left completely alone — reads cannot be added.
//          A platform is a population; a Target_Depth group is not. The 30x group
//          spans 18x (T7) to 35x (G400RS), so one group-wide fraction would
//          over-down-sample some platforms and under-down-sample others.
//
//   Stages (numbers are the publishDir prefixes under results/):
//     00  BUILD_MATRIX            define the refined_core sample x variant set
//                                 once, for every step that scores against it
//         DOWNSAMPLE_REGION       CRAM -> target-region BAMs, per-platform
//                                 fraction x 3 replicates (seeds 1/2/3).
//                                 Unpublished — see the process for why.
//     01  BUILD_FORCE_CALL_SITES  cohort VCF (03_annotate_af_norm_gt) ->
//                                 sites-only "alleles" VCF over the target regions
//     02  CALL_GVCF               HaplotypeCaller -ERC GVCF --alleles per BAM
//                                 (force-calls the cohort alleles into each GVCF)
//     03  JOINT_GENOTYPE          per rep: GenomicsDBImport + GenotypeGVCFs on
//                                 all samples -> ONE harmonized multi-sample VCF
//     04  NORMALIZE_VCF           split + left-align + set-id
//     05  COMPARE_TO_SITES        variant recovery vs the sites set + refined_core
//     06  PREP_CONCORDANCE_VCFS   per rep: the refined_core comparison matrix,
//                                 pre-QC + post-genotype-QC
//     07  GENOTYPE_CONCORDANCE    3 reps x pre/post vs the cohort's own genotypes
//                                 -> one log + figures
//     08  BUILD_ASSOC_GENOTYPES   THE DELIVERABLE: each model's cohort with the
//                                 down-sampled genotypes swapped in, ready to
//                                 re-run the association on
//     09  ASSOC_CONCORDANCE       what the swap changed, on the variants each
//                                 association actually reads
//
//   Joint genotyping (03) harmonizes alleles across all samples, so each variant
//   is a SINGLE record with no duplicates — unlike per-sample calling + bcftools
//   merge, which collapses differently-represented STR indels into duplicate IDs.
// ─────────────────────────────────────────────────────────────────────────────

// ── Inputs / references ──────────────────────────────────────────────────────
params.cram_info          = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tmp/cram.v6/cram.v6.summary.csv'
params.fasta              = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa'

// Per-platform keep fractions, measured from the CRAMs by tuning.fraction.nf.
// Columns: PLATFORM, TARGET_DEPTH, N, DEPTH_MEDIAN, READLEN, KEEP_FRACTION, NOTE.
// KEEP_FRACTION == 1 means "at or below baseline, do not touch".
params.platform_fractions = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/down_sampling/tuning.fraction/results/01_platform_fractions/platform_fractions.tsv'

// ── Outputs / runtime ────────────────────────────────────────────────────────
params.results_dir        = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/down_sampling/results'
params.script_dir         = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/down_sampling/scripts'
params.conda_env_activate = 'cteph_geno_pro'
params.samtools_threads   = 8

// ── Sample-sheet columns ─────────────────────────────────────────────────────
params.sample_id_col      = 'ID_JHRPv6'
params.platform_col       = 'WGS_Platform'
params.cram_path_col      = 'Cram_Path'
params.cram_found_col     = 'Cram_Found'

// ── Down-sampling spec ───────────────────────────────────────────────────────
// One or more target regions ("chr:start-end"), all merged into a single BAM per
// replicate. Override on the CLI with e.g.
//   --regions chr16:53703963-54121941,chr4:76306733-76311130
params.regions            = [
    'chr16:53703963-54121941',
    'chr4:76306733-76311130',
]
// One replicate per seed; the seed value is also the replicate index (rep1..3).
params.seeds              = [1, 2, 3]

// ── Force-calling spec ───────────────────────────────────────────────────────
// Per-chromosome cohort VCFs (analysis-ready, AF-annotated, genotype-normalized
// cteph_agp3k.v6 study cohort): ${cohort_vcf_dir}/<chr>.vqc.af.gtnorm.vcf.gz(.tbi).
// The alleles to force-call are the cohort FILTER=PASS variants within
// params.regions (sites-only).
params.cohort_vcf_dir     = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/03_annotate_af_norm_gt'
params.cohort_vcf_tpl     = 'CHR.vqc.af.gtnorm.vcf.gz'   // CHR is replaced by chr16, chr4, ...
params.gatk_sif           = '/home/b/b37974/simg/gatk_latest.sif'
params.gatk_java_opts     = '-Xmx8G -XX:ParallelGCThreads=4'
params.bcftools_threads   = 4

// ── The refined_core comparison matrix ───────────────────────────────────────
// The truth our down-sampled genotypes are compared against:
//
//   SAMPLES   refined_core_samples (FID<TAB>IID) — this list IS refined_core, and
//             is taken as given. Downstream steps match it BY ID against what they
//             actually hold and report what they scored.
//   VARIANTS  the target-region variants of refined_core_prefix.bim.
//   GENOTYPES refined_core_prefix (.bed/.bim/.fam).
//
// refined_core_prefix points at 14_fixed_model_prep/fixed_ready, NOT at
// 15_random_model_prep: fixed_ready keeps the RARE variants (3,612 in the target
// regions vs 976), which is exactly what this test needs — pull reads out and the
// calls that break first are the ones few reads supported.
//
// fixed_ready holds 2,166 of the 2,193 because the fixed-model prep prunes on
// relatedness (one of each related pair, so the fixed-effect model sees
// independent samples). Those 27 are still refined_core — they simply have no
// truth genotypes in this fileset. All 27 are HiSeqX 15x, the baseline platform,
// which is never down-sampled, so none of them could have been scored anyway.
params.refined_core_prefix  = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/14_fixed_model_prep/refined_core/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.fixed_ready'
params.refined_core_samples = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/PopGMM_output/refined_core.fid_iid.txt'

// ── The association genotypes ────────────────────────────────────────────────
// The deliverable: each base cohort with the down-sampled samples' re-called
// genotypes swapped in, so the association can be re-run on depth-levelled data.
// One fileset per model x replicate.
//
//   base      the fileset that model normally reads. It defines the sample and
//             variant universe; every sample that was not down-sampled keeps the
//             genotypes it already had.
//   region    the locus this model tests.
//   related   a list of related samples the model must NOT contain, or ''.
//             A fixed-effect model has no GRM to absorb relatedness, so its base
//             is pruned upstream; a random-effect model keeps them on purpose and
//             lets the GRM handle it. Declaring the list makes BUILD_ASSOC_
//             GENOTYPES verify the base rather than trust it — pointed at the
//             wrong fileset, a fixed-effect run would otherwise return a perfectly
//             normal-looking, wrong p-value.
params.wgs_results  = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results'
params.assoc_models = [
    [ name   : 'snp_based',
      base   : "${params.wgs_results}/15_random_model_prep/refined_core/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.random_model",
      region : 'chr16:53703963-54121941',
      related: '' ],
    [ name   : 'gene_based',
      base   : "${params.wgs_results}/14_fixed_model_prep/refined_core/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.maf_lt_threshold",
      region : 'chr4:76306733-76311130',
      related: "${params.wgs_results}/14_fixed_model_prep/refined_core/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.pihat_selected.exclude.fid_iid" ],
]

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

// The sample sheet has shipped both TAB- and COMMA-separated (the .csv name says
// nothing). Guessing wrong does not throw: every field lands in column 0, every
// lookup returns null, the filter matches nothing, and the run "succeeds" having
// down-sampled zero samples. So sniff the header and fail loudly.
def readTable(path, required) {
    def lines = file(path, checkIfExists: true).readLines().findAll { l -> l.trim() }
    if (!lines) error "Empty table: ${path}"
    def sep_re  = lines[0].count('\t') >= lines[0].count(',') && lines[0].contains('\t') ? /\t/ : /,/
    def hdr     = lines[0].split(sep_re, -1)*.trim()
    def missing = required.findAll { c -> !hdr.contains(c) }
    if (missing) {
        error "Table ${path} is missing column(s) ${missing}. Found: ${hdr}"
    }
    return lines.drop(1).collect { l ->
        [hdr, l.split(sep_re, -1)].transpose().collectEntries { k, v -> [(k): v?.trim()] }
    }
}

def regionList() {
    def rl = (params.regions instanceof List ? params.regions : "${params.regions}".split(','))
    return rl.collect { r -> r.toString().trim() }.findAll { r -> r }
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: extract the target regions once, then subsample to this sample's
//   platform fraction, three times.
// ─────────────────────────────────────────────────────────────────────────────
process DOWNSAMPLE_REGION {

    tag "${sample_id} | x${keep_fraction}"

    executor 'slurm'
    queue    'gr10478b'
    time     '12h'

    // No publishDir: the down-sampled BAMs are intermediates consumed by
    // CALL_GVCF through the channel. Publishing 6 glob-matched files per task at
    // high concurrency triggers a Nextflow 26.04 race
    // (ConcurrentModificationException in publishDir finalization) that aborts
    // the run, so we keep them in the work dir only.

    input:
    tuple val(sample_id), val(keep_fraction), path(cram), path(crai)

    output:
    tuple val(sample_id), path("${sample_id}.rep*.bam"), path("${sample_id}.rep*.bam.bai"), emit: bams

    script:
    def seed_list = params.seeds.join(' ')
    def region_args = regionList().join(' ')
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    # 1) Extract ALL target regions from the CRAM in a single pass.
    #    CRAM is coordinate-sorted, so the region output is already sorted.
    samtools view -@ ${params.samtools_threads} -b \\
        -T ${params.fasta} \\
        ${cram} \\
        ${region_args} \\
        -o region.full.bam
    samtools index region.full.bam

    # 2) Keep this platform's fraction of the reads, once per seed.
    #    --subsample/--subsample-seed rather than the older '-s INT.FRAC' spelling:
    #    there the fraction is the literal decimal string after the dot, so 0.05
    #    must be written '-s 1.050' and '-s 1.05' silently means 0.05 while
    #    '-s 1.5' means 0.5 — a fraction built by string concatenation is one
    #    typo away from a 10x error. The float form cannot be misread.
    #    Mates share a read name, so pairs are kept or dropped together.
    for seed in ${seed_list}; do
        samtools view -@ ${params.samtools_threads} -b \\
            --subsample ${keep_fraction} \\
            --subsample-seed \${seed} \\
            region.full.bam \\
            -o ${sample_id}.rep\${seed}.bam
        samtools index ${sample_id}.rep\${seed}.bam
    done
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: build the force-call "alleles" VCF from the cohort VCF
//   For every target region, pull the cohort variants (sites only) from the
//   matching per-chromosome VCF, then concat + coordinate-sort into one VCF.
// ─────────────────────────────────────────────────────────────────────────────
process BUILD_FORCE_CALL_SITES {

    tag "force_call_sites"

    executor 'slurm'
    queue    'gr10478b'
    time     '12h'

    publishDir "${params.results_dir}/01_force_call_sites", mode: 'symlink'

    input:
    val regions

    output:
    tuple path("force_call_sites.vcf.gz"), path("force_call_sites.vcf.gz.tbi"), emit: sites

    script:
    def region_list = (regions instanceof List ? regions : "${regions}".split(','))
    def region_args = region_list.collect { r -> r.toString().trim() }.findAll { r -> r }.join(' ')
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    : > sites_list.txt
    i=0
    for reg in ${region_args}; do
        chr=\${reg%%:*}
        cohort=${params.cohort_vcf_dir}/\$(echo ${params.cohort_vcf_tpl} | sed "s/CHR/\${chr}/")
        out=region_\${i}.sites.vcf.gz
        # sites-only (-G), FILTER=PASS (-f PASS), restricted to the region via .tbi
        bcftools view -G -f PASS -r \${reg} \${cohort} --threads ${params.bcftools_threads} -Oz -o \${out}
        bcftools index -t \${out}
        echo \${out} >> sites_list.txt
        i=\$((i+1))
    done

    # Merge per-region site sets in genomic order (-a: allow cross-file ordering),
    # then coordinate-sort to be safe against header contig-order quirks.
    bcftools concat -a -f sites_list.txt --threads ${params.bcftools_threads} -Oz -o sites.concat.vcf.gz
    bcftools sort sites.concat.vcf.gz -Oz -o force_call_sites.vcf.gz
    bcftools index -t force_call_sites.vcf.gz
    rm -f region_*.sites.vcf.gz* sites.concat.vcf.gz
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: per-sample GVCF calling on one down-sampled BAM.
//   HaplotypeCaller -ERC GVCF emits the reference-confidence model; --alleles
//   force-calls the cohort sites so their (cohort-harmonized) alleles are written
//   into every sample's GVCF. -L restricts to the target regions. Output filenames
//   keep the sample name from the BAM @RG SM tag; two explicit files per task
//   (g.vcf.gz + .tbi, not a glob) make publishing safe under the Nextflow 26.04
//   publishDir race.
// ─────────────────────────────────────────────────────────────────────────────
process CALL_GVCF {

    tag "${sample_id} | ${rep}"

    executor 'slurm'
    queue    'gr10478b'
    time     '12h'

    publishDir "${params.results_dir}/02_gvcf", mode: 'symlink'

    input:
    tuple val(sample_id), val(rep), path(bam), path(bai)
    tuple path(sites_vcf), path(sites_tbi)

    output:
    tuple val(rep),
          path("${sample_id}.${rep}.g.vcf.gz"),
          path("${sample_id}.${rep}.g.vcf.gz.tbi"), emit: gvcf

    script:
    def l_args = regionList().collect { r -> "-L ${r}" }.join(' ')
    def out_gvcf = "${sample_id}.${rep}.g.vcf.gz"
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    singularity exec \\
        --bind /LARGE0:/LARGE0 \\
        --bind /LARGE1:/LARGE1 \\
        ${params.gatk_sif} gatk --java-options "${params.gatk_java_opts}" HaplotypeCaller \\
        -R ${params.fasta} \\
        -I ${bam} \\
        ${l_args} \\
        -ERC GVCF \\
        --alleles ${sites_vcf} \\
        -O ${out_gvcf}

    # GATK already writes a .tbi; -f makes re-indexing idempotent on -resume.
    bcftools index -f -t ${out_gvcf}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: joint-genotype ALL samples of ONE replicate.
//   GenomicsDBImport harmonizes alleles across samples into a single workspace,
//   then GenotypeGVCFs emits ONE coordinate-sorted multi-sample VCF. Alleles are
//   consistent across samples, so no duplicate variant records arise (unlike
//   per-sample calling + bcftools merge). Restricted to the target regions.
//
//   Every cohort site is emitted regardless of evidence:
//     --force-output-intervals <sites>  force output at all cohort positions
//     -stand-call-conf 0                emission is decided by --alleles (baked
//                                       into the GVCFs), not by a QUAL threshold
//   so a cohort variant that all down-sampled samples lost still appears (all 0/0)
//   and its false-negatives are counted in the concordance.
// ─────────────────────────────────────────────────────────────────────────────
process JOINT_GENOTYPE {

    tag "${rep}"

    executor 'slurm'
    queue    'gr10478b'
    time     '12h'

    publishDir "${params.results_dir}/03_joint_vcf", mode: 'symlink'

    input:
    tuple val(rep), path(gvcfs), path(tbis)
    tuple path(sites_vcf), path(sites_tbi)

    output:
    tuple val(rep), path("${rep}.joint.vcf.gz"), path("${rep}.joint.vcf.gz.tbi"), emit: joint

    script:
    def l_args = regionList().collect { r -> "-L ${r}" }.join(' ')
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    # sample-name-map: internal sample name (from the GVCF) -> GVCF path
    : > map.txt
    for g in *.g.vcf.gz; do
        s=\$(bcftools query -l "\$g")
        printf '%s\\t%s\\n' "\$s" "\$g" >> map.txt
    done

    # 1) Import all per-sample GVCFs into one GenomicsDB workspace (allele-harmonized)
    rm -rf gdb
    singularity exec \\
        --bind /LARGE0:/LARGE0 \\
        --bind /LARGE1:/LARGE1 \\
        ${params.gatk_sif} gatk --java-options "${params.gatk_java_opts}" GenomicsDBImport \\
        --sample-name-map map.txt \\
        --genomicsdb-workspace-path gdb \\
        --batch-size 50 \\
        ${l_args}

    # 2) Joint-genotype -> one harmonized, duplicate-free multi-sample VCF
    singularity exec \\
        --bind /LARGE0:/LARGE0 \\
        --bind /LARGE1:/LARGE1 \\
        ${params.gatk_sif} gatk --java-options "${params.gatk_java_opts}" GenotypeGVCFs \\
        -R ${params.fasta} \\
        -V gendb://gdb \\
        ${l_args} \\
        --force-output-intervals ${sites_vcf} \\
        -stand-call-conf 0 \\
        -O ${rep}.joint.vcf.gz

    bcftools index -f -t ${rep}.joint.vcf.gz
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: PREPARE_VCF-style representation normalization (necessary parts only).
//   Mirrors select.auto.par.v6.nf::PREPARE_VCF for the steps that matter here:
//     - split multiallelics + left-align + normalize (--check-ref s fixes REF)
//     - drop non-variant placeholder ALTs (* and .)   [see below]
//     - set ID = CHROM:POS:REF:ALT
//   The sample-selection / -f PASS / MAC>=1 steps are intentionally omitted:
//   inputs are already our samples, joint calls are unfiltered (FILTER='.'), and
//   we keep every variant to preserve the maximal count.
//   Applied to each joint replicate VCF AND to the sites VCF so both share one key.
//
//   Why ALT='*' and ALT='.' appear and are removed (both are NON-variants):
//     ALT='*'  spanning-deletion placeholder: an upstream deletion overlaps this
//              position, so GATK marks it with '*'. Not a variant at this site.
//     ALT='.'  GenotypeGVCFs --force-output-intervals forces a record at every
//              cohort position; where a replicate has no genotypable ALT (no/low
//              coverage, or the forced --alleles allele isn't supported) it emits
//              a reference-only row with ALT='.'. Not a variant.
//   Neither carries a real CHROM:POS:REF:ALT, so they only inflate "extra" counts.
// ─────────────────────────────────────────────────────────────────────────────
process NORMALIZE_VCF {

    tag "${label}"

    executor 'slurm'
    queue    'gr10478b'
    time     '12h'

    publishDir "${params.results_dir}/04_norm_vcf", mode: 'symlink'

    input:
    tuple val(label), path(vcf), path(tbi)

    output:
    tuple val(label), path("${label}.norm.vcf.gz"), path("${label}.norm.vcf.gz.tbi"), emit: norm

    script:
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    bcftools norm \\
        --multiallelics -any \\
        --fasta-ref ${params.fasta} \\
        --check-ref s \\
        --threads ${params.bcftools_threads} ${vcf} -Ou | \\
      bcftools filter --threads ${params.bcftools_threads} -e 'ALT="*" || ALT="."' -Ou | \\
      bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' \\
        --threads ${params.bcftools_threads} -Oz -o ${label}.norm.vcf.gz
    bcftools index -t ${label}.norm.vcf.gz
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: compare a normalized replicate VCF against the normalized sites VCF
//   and write a log (variant counts + recovered / missed / extra breakdown).
//   Delegates the set comparison to scripts/compare_to_sites.py.
// ─────────────────────────────────────────────────────────────────────────────
process COMPARE_TO_SITES {

    tag "${rep}"

    executor 'slurm'
    queue    'gr10478b'
    time     '4h'

    publishDir "${params.results_dir}/05_compare_log", mode: 'symlink'

    input:
    tuple val(rep), path(rep_vcf), path(rep_tbi)
    tuple path(sites_vcf), path(sites_tbi)
    path refined_ids
    path script

    output:
    path "${rep}.compare.log", emit: log
    path "${rep}.*.txt",       emit: diffs

    script:
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    python3 ${script} \\
        --rep-vcf ${rep_vcf} \\
        --sites-vcf ${sites_vcf} \\
        --refined-ids ${refined_ids} \\
        --label ${rep} \\
        --out-log ${rep}.compare.log \\
        --out-dir .
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: build the refined_core comparison matrix for ONE replicate, pre- and
//   post-genotype-QC.
//     pre_qc  = replicate genotypes on the matrix (samples x variants), BEFORE
//               genotype QC — this is what the down-sampled calls look like raw.
//     post_qc = the same after adding FORMAT/AF (AlleleFraction), unphasing, and
//               FILTER_GENOTYPE-style genotype QC (GQ<20 / DP<8 / allele-balance
//               -> no-call). Mirrors PREPARE_VCF's downstream
//               ANNOTATE_AF_NORM_GT + FILTER_GENOTYPE.
//   Both VCFs are cut to the matrix here, so the published artefact IS the
//   comparison set rather than something the next step has to re-derive.
// ─────────────────────────────────────────────────────────────────────────────
process PREP_CONCORDANCE_VCFS {

    tag "${rep}"

    executor 'slurm'
    queue    'gr10478b'
    time     '6h'

    publishDir "${params.results_dir}/06_concordance_vcf", mode: 'symlink'

    input:
    tuple val(rep), path(rep_vcf), path(rep_tbi)
    path refined_ids
    path matrix_samples
    path script
    path backfill_script

    output:
    tuple val(rep),
          path("${rep}.pre_qc.vcf.gz"),  path("${rep}.pre_qc.vcf.gz.tbi"),
          path("${rep}.post_qc.vcf.gz"), path("${rep}.post_qc.vcf.gz.tbi"), emit: vcfs

    script:
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    bash ${script} \\
        ${rep} \\
        ${rep_vcf} \\
        ${refined_ids} \\
        ${matrix_samples} \\
        ${backfill_script} \\
        ${params.fasta} \\
        ${params.gatk_sif} \\
        "${params.gatk_java_opts}"
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: genotype concordance of ALL replicates (pre_qc & post_qc) vs the
//   cohort's own genotypes, on the refined_core matrix. One combined log + two
//   figures (3 reps x pre/post confusion matrices with concordance / FPR / FNR /
//   miss-rate, plus smiss/vmiss). Truth is exported as ALT dosage
//   (--export-allele forces ALT); test GT come straight from the VCFs.
// ─────────────────────────────────────────────────────────────────────────────
process GENOTYPE_CONCORDANCE {

    tag "all_reps"

    executor 'slurm'
    queue    'gr10478b'
    time     '6h'

    publishDir "${params.results_dir}/07_genotype_concordance", mode: 'symlink'

    input:
    path vcfs
    path refined_ids
    path matrix_samples
    path sh_script
    path py_script

    output:
    path "genotype_concordance.log",            emit: log
    path "genotype_concordance_matrices.png",   emit: figure
    path "genotype_missingness.png",            emit: figure_missing

    script:
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    bash ${sh_script} \\
        ${params.refined_core_prefix} \\
        ${refined_ids} \\
        ${matrix_samples} \\
        ${py_script}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: define the refined_core comparison matrix once, up front.
//   Emits the variant IDs and the sample IDs that every downstream step shares,
//   plus a log of how the intersection was reached. Doing this once — rather
//   than re-deriving it inside three different scripts from the same .bim — means
//   the matrix cannot quietly differ between steps.
// ─────────────────────────────────────────────────────────────────────────────
process BUILD_MATRIX {

    tag "refined_core"

    executor 'slurm'
    queue    'gr10478b'
    time     '2h'

    publishDir "${params.results_dir}/00_matrix", mode: 'symlink'

    input:
    path core_samples
    path script

    output:
    path "matrix.variants.txt", emit: variants
    path "matrix.samples.txt",  emit: samples
    path "matrix.log",          emit: log

    script:
    def region_args = regionList().join(' ')
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    bash ${script} \\
        ${params.refined_core_prefix} \\
        ${core_samples} \\
        matrix.variants.txt \\
        matrix.samples.txt \\
        matrix.log \\
        ${region_args}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: the deliverable — one model's base cohort with ONE replicate's
//   re-called genotypes swapped in for the down-sampled samples. Run the
//   association on this: the depth difference it was accused of exploiting is
//   gone, so an effect that survives cannot be made of it.
// ─────────────────────────────────────────────────────────────────────────────
process BUILD_ASSOC_GENOTYPES {

    tag "${rep} | ${model.name}"

    executor 'slurm'
    queue    'gr10478b'
    time     '4h'

    publishDir "${params.results_dir}/08_assoc_genotypes", mode: 'symlink'

    input:
    tuple val(rep), path(post_qc), path(post_tbi), val(model)
    path script

    output:
    tuple val(model.name), val(rep),
          path("${rep}.${model.name}.bed"),
          path("${rep}.${model.name}.bim"),
          path("${rep}.${model.name}.fam"),
          path("${rep}.${model.name}.ds_samples.txt"), emit: geno
    path "${rep}.${model.name}.log", emit: log

    script:
    def related = model.related ?: 'NONE'
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    bash ${script} \\
        ${rep} \\
        ${model.name} \\
        ${model.base} \\
        ${model.region} \\
        ${post_qc} \\
        ${rep}.${model.name} \\
        ${related}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: score every association fileset against its base, on the variants the
//   association actually reads and the samples that actually changed. 07 answers
//   "did re-calling reproduce the cohort" on the 3,612-variant matrix; this
//   answers it on the ~970 and ~22 variants the two tests really use, which is
//   where a moved odds ratio would have to come from.
// ─────────────────────────────────────────────────────────────────────────────
process ASSOC_CONCORDANCE {

    tag "all_models"

    executor 'slurm'
    queue    'gr10478b'
    time     '4h'

    publishDir "${params.results_dir}/09_assoc_concordance", mode: 'symlink'

    input:
    path filesets
    path spec
    path sh_script
    path py_script
    // Staged so `from genotype_concordance import ...` resolves: the metrics are
    // imported, not restated, so the two reports cannot drift apart.
    path cc_py_script

    output:
    path "assoc_concordance.log", emit: log
    path "assoc_concordance.png", emit: figure

    script:
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    bash ${sh_script} ${py_script} ${spec}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Workflow
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // Scripts as path inputs, not bare param paths: only a path input is
    // content-hashed, so editing one actually invalidates its task on -resume.
    matrix_script  = file("${params.script_dir}/build_matrix.sh",           checkIfExists: true)
    compare_script = file("${params.script_dir}/compare_to_sites.py",       checkIfExists: true)
    prep_script    = file("${params.script_dir}/prep_concordance_vcfs.sh",  checkIfExists: true)
    backfill_script = file("${params.script_dir}/backfill_rgq_to_gq.py",    checkIfExists: true)
    cc_sh_script   = file("${params.script_dir}/genotype_concordance.sh",   checkIfExists: true)
    cc_py_script   = file("${params.script_dir}/genotype_concordance.py",   checkIfExists: true)
    assoc_script       = file("${params.script_dir}/build_assoc_genotypes.sh", checkIfExists: true)
    assoc_cc_sh_script = file("${params.script_dir}/assoc_concordance.sh",     checkIfExists: true)
    assoc_cc_py_script = file("${params.script_dir}/assoc_concordance.py",     checkIfExists: true)

    // ── Sample selection: platforms with KEEP_FRACTION < 1 ────────────────────
    def frac_rows = readTable(params.platform_fractions, ['PLATFORM', 'KEEP_FRACTION'])
    def frac_of   = frac_rows.collectEntries { r -> [(r['PLATFORM']): (r['KEEP_FRACTION'] as double)] }
    def to_sample = frac_of.findAll { _pf, f -> f < 1.0 }

    log.info "[down_sampling] per-platform keep fractions (${params.platform_fractions}):"
    frac_of.sort { e -> -e.value }.each { pf, f ->
        log.info "[down_sampling]   ${pf.padRight(22)} ${String.format('%.4f', f)}" +
                 (f < 1.0 ? "  -> down-sample" : "  -> untouched (at or below baseline)")
    }
    if (!to_sample) {
        error "Every platform in ${params.platform_fractions} has KEEP_FRACTION == 1, so there " +
              "is nothing to down-sample. Re-run tuning.fraction.nf, or check its baseline."
    }

    def rows = readTable(params.cram_info,
        [params.sample_id_col, params.platform_col, params.cram_path_col, params.cram_found_col])

    // A platform present in the sheet but absent from the fractions table has no
    // fraction, so it would be dropped without anyone noticing. Refuse instead.
    def sheet_platforms = rows.collect { r -> r[params.platform_col] }.unique().findAll { p -> p }
    def unknown = sheet_platforms.findAll { p -> !frac_of.containsKey(p) }
    if (unknown) {
        error "Platform(s) ${unknown} appear in ${params.cram_info} but not in " +
              "${params.platform_fractions}, so they have no keep fraction. Re-run " +
              "tuning.fraction.nf against this sheet before down-sampling."
    }

    def selected = rows.findAll { r ->
        to_sample.containsKey(r[params.platform_col]) && r[params.cram_found_col] == 'True'
    }
    if (!selected) {
        error "No sample of a down-sampled platform has a CRAM in ${params.cram_info} — " +
              "refusing to run a pipeline that would down-sample nothing."
    }
    log.info "[down_sampling] ${selected.size()} samples to down-sample:"
    selected.groupBy { r -> r[params.platform_col] }.sort { e -> -e.value.size() }.each { pf, rs ->
        log.info "[down_sampling]   ${pf.padRight(22)} ${rs.size().toString().padLeft(4)} samples " +
                 "x ${String.format('%.4f', to_sample[pf])}"
    }

    cram_ch = channel
        .fromList(selected)
        .map { row ->
            tuple(row[params.sample_id_col],
                  to_sample[row[params.platform_col]],
                  file(row[params.cram_path_col]),
                  file("${row[params.cram_path_col]}.crai"))
        }

    // ── The comparison matrix, defined once and shared ────────────────────────
    matrix = BUILD_MATRIX(file(params.refined_core_samples, checkIfExists: true), matrix_script)

    // ── 00: down-sample ───────────────────────────────────────────────────────
    ds_out = DOWNSAMPLE_REGION(cram_ch)

    // ── 01: the shared force-call alleles VCF (one job for all samples) ───────
    sites_ch = BUILD_FORCE_CALL_SITES(channel.value(params.regions)).sites

    // Flatten (sample, [bams], [bais]) -> one (sample, rep, bam, bai) per replicate.
    per_rep_ch = ds_out.bams.flatMap { sid, bams, bais ->
        def blist = (bams instanceof List ? bams : [bams]).sort { a, b -> a.name <=> b.name }
        def ilist = (bais instanceof List ? bais : [bais]).sort { a, b -> a.name <=> b.name }
        [blist, ilist].transpose().collect { pair ->
            def bam = pair[0]
            def bai = pair[1]
            // The name is "<sample>.<rep>.bam", so the label is the token before .bam.
            def parts = bam.name.tokenize('.')
            def rep   = parts.size() >= 3 ? parts[-2] : ''
            if (!(rep ==~ /rep\d+/)) {
                // Falling back to the basename would silently invent a label like
                // 'PHOM0001.rep1' and split each joint call into one-sample groups.
                // Better to stop than to emit a plausible wrong answer.
                error "Cannot read the replicate label from BAM '${bam.name}'. " +
                      "Expected <sample>.rep<N>.bam."
            }
            tuple(sid, rep, bam, bai)
        }
    }

    // ── 02: per-sample GVCF calling (cohort alleles forced in) ────────────────
    // sites_ch is a value channel (BUILD_FORCE_CALL_SITES was fed a value input),
    // so it is reused across all per-replicate items automatically.
    gvcf_out = CALL_GVCF(per_rep_ch, sites_ch)

    // ── 03: joint-genotype per replicate ─────────────────────────────────────
    by_rep_ch = gvcf_out.gvcf.groupTuple()               // (rep, [gvcf...], [tbi...])
    joint_ch  = JOINT_GENOTYPE(by_rep_ch, sites_ch).joint // (rep, joint_vcf, tbi)

    // ── 04: normalize the replicate VCFs AND the sites VCF the same way ──────
    norm_in_ch = joint_ch
        .mix( sites_ch.map { vcf, tbi -> tuple('sites', vcf, tbi) } )
    norm_ch = NORMALIZE_VCF(norm_in_ch).norm  // (label, norm_vcf, tbi)

    sites_norm_ch = norm_ch
        .filter { label, _v, _t -> label == 'sites' }
        .map    { _label, v, t -> tuple(v, t) }
    reps_norm_ch = norm_ch
        .filter { label, _v, _t -> label != 'sites' }

    // ── 05: variant recovery per replicate ───────────────────────────────────
    COMPARE_TO_SITES(reps_norm_ch, sites_norm_ch.first(), matrix.variants, compare_script)

    // ── 06/07: the matrix pre/post genotype QC, then one combined concordance ─
    prep_ch = PREP_CONCORDANCE_VCFS(reps_norm_ch, matrix.variants, matrix.samples,
                                    prep_script, backfill_script).vcfs
    cc_inputs_ch = prep_ch
        .flatMap { _rep, pre, pre_tbi, post, post_tbi -> [pre, pre_tbi, post, post_tbi] }
        .collect()
    GENOTYPE_CONCORDANCE(cc_inputs_ch, matrix.variants, matrix.samples,
                         cc_sh_script, cc_py_script)

    // ── 08: the association genotypes, one per model x replicate ─────────────
    assoc_in_ch = prep_ch
        .map { rep, _pre, _pre_tbi, post, post_tbi -> tuple(rep, post, post_tbi) }
        .combine( channel.fromList(params.assoc_models) )
    assoc_ch = BUILD_ASSOC_GENOTYPES(assoc_in_ch, assoc_script).geno

    // ── 09: what the swap actually changed, where the association looks ──────
    // The spec pairs each fileset with the base it must be scored against; the
    // script cannot infer that from the staged files alone.
    assoc_spec_ch = assoc_ch
        .map { model, rep, _bed, _bim, _fam, _ds ->
            def base = params.assoc_models.find { m -> m.name == model }.base
            "${model}\t${rep}\t${base}\n"
        }
        .collectFile(name: 'assoc_spec.tsv', sort: true)

    assoc_files_ch = assoc_ch
        .flatMap { _model, _rep, bed, bim, fam, ds -> [bed, bim, fam, ds] }
        .collect()

    ASSOC_CONCORDANCE(assoc_files_ch, assoc_spec_ch,
                      assoc_cc_sh_script, assoc_cc_py_script, cc_py_script)
}
