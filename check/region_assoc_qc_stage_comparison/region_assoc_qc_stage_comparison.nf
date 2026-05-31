nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// Parameters: Genomic Region of Interest
// ─────────────────────────────────────────────────────────────────────────────
params.genomic_locus                = "chr16:53703963-54121941"
params.target_loci                  = "chr16:53887925:T:C"

// ─────────────────────────────────────────────────────────────────────────────
// Parameters: Input Data (QC Stage Sources)
// ─────────────────────────────────────────────────────────────────────────────
params.pre_genotype_qc_vcf_dir      = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/02_filter_vqc"
params.post_genotype_qc_plink_dir   = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/05_vcf_to_plink"
params.post_variant_qc_plink_dir    = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/14_fixed_model_prep"
params.post_variant_qc_plink_prefix = "${params.post_variant_qc_plink_dir}/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.fixed_ready"

//-----------------------------------------------------------------------------
// PopGMM Configuration
//------------------------------------------------------------------------------
params.popgmm = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/PopGMM_output/mainland_subcluster_samples.fid_iid.txt'

// ─────────────────────────────────────────────────────────────────────────────
// Parameters: Association Analysis
// ─────────────────────────────────────────────────────────────────────────────
params.assoc_input_root             = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results"
params.pheno_rel_path               = "16_cov_pheno_prep/popgmm_subset_on_bbj_pcs.pheno.tsv"
params.pheno_name                   = "PHENO1"
params.covar_names                  = "SEX,PC1_AVG-PC10_AVG"
params.covar_pc_source              = "bbj"   // own | bbj
params.include_age_z                = params.covar_names
    .split(',')
    .collect { String s -> s.trim().toUpperCase() }
    .contains('AGE_Z')
params.covar_rel_path               = params.covar_pc_source == 'bbj' \
    ? (params.include_age_z \
        ? '16_cov_pheno_prep/popgmm_subset_on_bbj_pcs.cov.sex_age_agez.tsv' \
        : '16_cov_pheno_prep/popgmm_subset_on_bbj_pcs.cov.sex.tsv') \
    : (params.include_age_z \
        ? '16_cov_pheno_prep/popgmm_relatedness_aware_projection.cov.sex_age_agez.tsv' \
        : '16_cov_pheno_prep/popgmm_relatedness_aware_projection.cov.sex.tsv')
params.age_na_remove_rel_path       = "16_cov_pheno_prep/popgmm_subset_on_bbj_pcs.age_na.fid_iid"
params.firth_mode                   = "no-firth"   // no-firth | firth-fallback | firth
params.maf_filter_min               = 0.01          // minimum MAF; variants below this are excluded before association
params.maf_filter_source            = "ctrl"        // ctrl (MAF computed in controls only) | all (all samples)
params.output_prefix                = "cteph_agp3k.v6"

// ─────────────────────────────────────────────────────────────────────────────
// Parameters: Runtime
// ─────────────────────────────────────────────────────────────────────────────
params.bcftools_threads             = 4
params.plink2_threads               = 8
params.only_validate                = false
params.script_dir                   = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/region_assoc_qc_stage_comparison/scripts"
params.results_dir                  = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/region_assoc_qc_stage_comparison/results"
params.conda_env_activate           = "cteph_geno_pro"
params.gene_gtf                     = "/LARGE1/gr10478/platform/JHRPv6/workspace/pipeline/data/GTF/Homo_sapiens.GRCh38.97.chr.gtf"   // optional hg38 GTF(.gz) path for bottom gene-structure panel


// ─────────────────────────────────────────────────────────────────────────────
// Helper: parse "chrN:START-END" into a map
// ─────────────────────────────────────────────────────────────────────────────
def parseGenomicLocus(String locus) {
    def (chrom, range) = locus.tokenize(':')
    def (start, end)   = range.tokenize('-')
    return [chrom: chrom.trim(), start: start.trim(), end: end.trim()]
}

def validateParams() {
    if (!(params.genomic_locus ==~ /^chr?[0-9XYM]+:[0-9]+-[0-9]+$/)) {
        error "Invalid params.genomic_locus='${params.genomic_locus}'. Expected format like chr16:53703963-54121941"
    }

    def loc = parseGenomicLocus(params.genomic_locus)
    if (loc.start as long > loc.end as long) {
        error "Invalid params.genomic_locus='${params.genomic_locus}'. START must be <= END"
    }

    def validCovarPcSource = ['own', 'bbj']
    if (!validCovarPcSource.contains(params.covar_pc_source)) {
        error "Invalid params.covar_pc_source='${params.covar_pc_source}'. Allowed values: ${validCovarPcSource.join(', ')}"
    }

    def validFirthMode = ['no-firth', 'firth-fallback', 'firth']
    if (!validFirthMode.contains(params.firth_mode)) {
        error "Invalid params.firth_mode='${params.firth_mode}'. Allowed values: ${validFirthMode.join(', ')}"
    }

    if ((params.bcftools_threads as int) <= 0) {
        error "Invalid params.bcftools_threads='${params.bcftools_threads}'. Must be a positive integer"
    }

    if ((params.plink2_threads as int) <= 0) {
        error "Invalid params.plink2_threads='${params.plink2_threads}'. Must be a positive integer"
    }

    if (!(params.only_validate in [true, false])) {
        error "Invalid params.only_validate='${params.only_validate}'. Must be true or false"
    }

    def mafMin = params.maf_filter_min as double
    if (mafMin < 0 || mafMin >= 0.5) {
        error "Invalid params.maf_filter_min='${params.maf_filter_min}'. Must be in [0, 0.5)"
    }

    def validMafSource = ['ctrl', 'all']
    if (!validMafSource.contains(params.maf_filter_source)) {
        error "Invalid params.maf_filter_source='${params.maf_filter_source}'. Allowed values: ${validMafSource.join(', ')}"
    }

    def covarTagScript = file("${params.script_dir}/covar_tag_from_names.py")
    if (!covarTagScript.exists()) {
        error "Missing covariate tag script: ${covarTagScript}"
    }

    if (!params.popgmm) {
        error "Missing required params.popgmm. All QC stages must be subset using PopGMM intersection."
    }

    def popgmmFile = file(params.popgmm)
    if (!popgmmFile.exists()) {
        error "Missing PopGMM file: ${params.popgmm}. All QC stages must be subset using PopGMM intersection."
    }
}

def chromAliases(String chrom) {
    def c = chrom.trim()
    def num = c.replaceFirst(/^chr/, '')
    return ["chr${num}", num].unique()
}

def firstExistingPath(List<String> candidates, String label) {
    def found = candidates.find { p -> file(p).exists() }
    if (!found) {
        error "Unable to resolve ${label}. Tried: ${candidates.join(', ')}"
    }
    return file(found)
}

def requiredPath(String p, String label) {
    def fp = file(p)
    if (!fp.exists()) {
        error "Missing ${label}: ${p}"
    }
    return fp
}


// ─────────────────────────────────────────────────────────────────────────────
// Process 1: Pre-genotype-QC stage  (VCF → PLINK)
//   Source : {chrom}.vqc.vcf.gz  (02_filter_vqc)
//   Actions: subset genomic locus  →  retain reference sample set
//            →  remove monomorphic variants  →  PLINK2 BED/BIM/FAM
// ─────────────────────────────────────────────────────────────────────────────
process EXTRACT_REGION_PRE_GENOTYPE_QC {

    tag "pre_genotype_qc | ${params.genomic_locus}"

    executor 'slurm'
    queue    'gr10478b'
    time     '36h'

    publishDir "${params.results_dir}/00_pre_genotype_qc_region", mode: 'symlink'

    input:
    path vcf
    path vcf_tbi
    path ref_fam

    output:
    tuple path("pre_genotype_qc_region.bed"), path("pre_genotype_qc_region.bim"), path("pre_genotype_qc_region.fam"), emit: plink
    path "pre_genotype_qc_region.log", emit: log

    script:
    def loc   = parseGenomicLocus(params.genomic_locus)
    def chrom = loc.chrom
    def start = loc.start
    def end   = loc.end
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}

    # Extract IIDs from reference fam for bcftools sample subsetting
    awk '{print \$2}' ${ref_fam} > samples_iid.txt

    # Subset VCF to target locus and reference sample set
    bcftools view \\
        -r ${chrom}:${start}-${end} \\
        -S samples_iid.txt \\
        --threads ${params.bcftools_threads} \\
        ${vcf} \\
        -O z -o region_subset.vcf.gz

    tabix -p vcf region_subset.vcf.gz

    # Convert to PLINK2 BED; --mac 1 removes monomorphic variants
    # --double-id: sets FID = IID = VCF sample name
    plink2 \\
        --vcf region_subset.vcf.gz \\
        --double-id \\
        --mac 1 \\
        --make-bed \\
        --threads ${params.plink2_threads} \\
        --out pre_genotype_qc_region
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Process 2: Post-genotype-QC stage  (PLINK → PLINK)
//   Source : {chrom}.plink.bed/bim/fam  (05_vcf_to_plink)
//   Actions: subset genomic locus  →  retain reference sample set (--keep)
//            →  remove monomorphic variants  →  PLINK2 BED/BIM/FAM
// ─────────────────────────────────────────────────────────────────────────────
process EXTRACT_REGION_POST_GENOTYPE_QC {

    tag "post_genotype_qc | ${params.genomic_locus}"

    executor 'slurm'
    queue    'gr10478b'
    time     '36h'

    publishDir "${params.results_dir}/00_post_genotype_qc_region", mode: 'symlink'

    input:
    path bed
    path bim
    path fam
    path ref_fam

    output:
    tuple path("post_genotype_qc_region.bed"), path("post_genotype_qc_region.bim"), path("post_genotype_qc_region.fam"), emit: plink
    path "post_genotype_qc_region.log", emit: log

    script:
    def loc       = parseGenomicLocus(params.genomic_locus)
    def chrom     = loc.chrom
    def chrom_num = chrom.replaceFirst(/^chr/, '')
    def start     = loc.start
    def end       = loc.end
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}

    # Auto-detect whether the BIM uses "chr*" or "*" chromosome naming
    bim_chrom=\$(awk 'NR==1{print \$1}' ${bim})
    if [[ "\${bim_chrom}" == chr* ]]; then
        chr_arg="${chrom}"
    else
        chr_arg="${chrom_num}"
    fi

    plink2 \\
        --bfile ${bed.baseName} \\
        --chr \${chr_arg} \\
        --from-bp ${start} \\
        --to-bp ${end} \\
        --keep ${ref_fam} \\
        --mac 1 \\
        --make-bed \\
        --threads ${params.plink2_threads} \\
        --out post_genotype_qc_region
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Process 3: Post-variant-QC stage  (PLINK → PLINK)
//   Source : post_variant_qc_plink_prefix.bed/bim/fam  (14_fixed_model_prep)
//   Actions: subset genomic locus  →  remove monomorphic variants
//            →  PLINK2 BED/BIM/FAM
//   Note   : sample set already finalised in the input fam; no --keep required
// ─────────────────────────────────────────────────────────────────────────────
process EXTRACT_REGION_POST_VARIANT_QC {

    tag "post_variant_qc | ${params.genomic_locus}"

    executor 'slurm'
    queue    'gr10478b'
    time     '36h'

    publishDir "${params.results_dir}/00_post_variant_qc_region", mode: 'symlink'

    input:
    path bed
    path bim
    path fam

    output:
    tuple path("post_variant_qc_region.bed"), path("post_variant_qc_region.bim"), path("post_variant_qc_region.fam"), emit: plink
    path "post_variant_qc_region.log", emit: log

    script:
    def loc       = parseGenomicLocus(params.genomic_locus)
    def chrom     = loc.chrom
    def chrom_num = chrom.replaceFirst(/^chr/, '')
    def start     = loc.start
    def end       = loc.end
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}

    # Auto-detect whether the BIM uses "chr*" or "*" chromosome naming
    bim_chrom=\$(awk 'NR==1{print \$1}' ${bim})
    if [[ "\${bim_chrom}" == chr* ]]; then
        chr_arg="${chrom}"
    else
        chr_arg="${chrom_num}"
    fi

    plink2 \\
        --bfile ${bed.baseName} \\
        --chr \${chr_arg} \\
        --from-bp ${start} \\
        --to-bp ${end} \\
        --mac 1 \\
        --make-bed \\
        --threads ${params.plink2_threads} \\
        --out post_variant_qc_region
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Process 4: Optional PopGMM intersection subset (per QC stage)
//   keep samples = intersection(ref_fam, params.popgmm)
//   Applied after region extraction and before association/LD.
// ─────────────────────────────────────────────────────────────────────────────
process SUBSET_REGION_BY_POPGMM {

    tag "${qc_stage} | popgmm_intersection"

    executor 'slurm'
    queue    'gr10478b'
    time     '36h'

    publishDir "${params.results_dir}/00_popgmm_intersection_region", mode: 'symlink'

    input:
    tuple val(qc_stage), path(bed), path(bim), path(fam), path(ref_fam), path(popgmm)

    output:
    tuple val(qc_stage), path("subset_region.bed"), path("subset_region.bim"), path("subset_region.fam"), emit: plink
    path "subset_region.log", emit: log

    script:
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}

    # Build keep list from shared (FID, IID) between popgmm and ref_fam.
    awk 'BEGIN{OFS="\\t"}
         FNR==NR {
             if (NF < 2) next
             if (toupper(\$1)=="FID" && toupper(\$2)=="IID") next
             key[\$1 OFS \$2]=1
             next
         }
         {
             if ((\$1 OFS \$2) in key) print \$1, \$2
         }' ${popgmm} ${ref_fam} > keep_intersection.txt

    if [[ ! -s keep_intersection.txt ]]; then
        echo "Error: PopGMM/ref_fam intersection is empty for stage ${qc_stage}" >&2
        exit 1
    fi

    echo "qc_stage=${qc_stage}" > subset_region.log
    echo "n_keep=\$(wc -l < keep_intersection.txt)" >> subset_region.log

    plink2 \\
        --bfile ${bed.baseName} \\
        --keep keep_intersection.txt \\
        --mac 1 \\
        --make-bed \\
        --threads ${params.plink2_threads} \\
        --out subset_region
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Process 4: Build shared association variant list from pre-genotype-QC only
//   Strategy:
//   - Compute frequency from pre_genotype_qc_region
//   - Apply MAF rule once (params.maf_filter_source: ctrl | all)
//   - Reuse the resulting variant ID list for all three association runs
// ─────────────────────────────────────────────────────────────────────────────
process BUILD_ASSOC_VARIANT_LIST {

    tag "variant_list | ${params.genomic_locus} | ${params.maf_filter_source}"

    executor 'slurm'
    queue    'gr10478b'
    time     '36h'

    publishDir "${params.results_dir}/01_region_assoc/variant_list", mode: 'symlink'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    path "assoc_variant_ids.txt", emit: variant_ids
    path "variant_list_build.log", emit: log

    script:
    def pheno_file = "${params.assoc_input_root}/${params.pheno_rel_path}"
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}

    # Keep a concise execution trace for debugging
    {
        echo "maf_filter_source=${params.maf_filter_source}"
        echo "maf_filter_min=${params.maf_filter_min}"
        echo "bfile=${bed.baseName}"
    } > variant_list_build.log

    if [[ "${params.maf_filter_source}" == "ctrl" ]]; then
        awk 'NR==1 {
                 for (i=1; i<=NF; i++) if (\$i == "${params.pheno_name}") col=i
                 if (!col) {
                     print "Error: phenotype column ${params.pheno_name} not found in " FILENAME > "/dev/stderr"
                     exit 2
                 }
                 next
             }
             \$col == 1 { print \$1, \$2 }' ${pheno_file} > ctrl_samples.txt

        if [[ ! -s ctrl_samples.txt ]]; then
            echo "Error: no control samples found in ${pheno_file} for ${params.pheno_name}==1" >&2
            exit 1
        fi

        plink2 \
            --bfile ${bed.baseName} \
            --keep ctrl_samples.txt \
            --freq \
            --threads ${params.plink2_threads} \
            --out maf_source_freq
    else
        plink2 \
            --bfile ${bed.baseName} \
            --freq \
            --threads ${params.plink2_threads} \
            --out maf_source_freq
    fi

    awk 'NR==1 {
             for (i=1; i<=NF; i++) h[\$i]=i
             if (!("ID" in h) || !("ALT_FREQS" in h)) {
                 print "Error: expected ID and ALT_FREQS columns in " FILENAME > "/dev/stderr"
                 exit 2
             }
             id_col = h["ID"]
             af_col = h["ALT_FREQS"]
             next
         }
         {
             n = split(\$af_col, a, ",")
             maf_min = 1.0
             for (j=1; j<=n; j++) {
                 f = a[j] + 0
                 maf = (f <= 0.5) ? f : 1 - f
                 if (maf < maf_min) maf_min = maf
             }
             if (maf_min >= ${params.maf_filter_min}) print \$id_col
         }' maf_source_freq.afreq > assoc_variant_ids.txt

    if [[ ! -s assoc_variant_ids.txt ]]; then
        echo "Error: shared MAF filter kept 0 variants (source=${params.maf_filter_source}, min=${params.maf_filter_min})." >&2
        exit 1
    fi

    echo "n_variant_ids=\$(wc -l < assoc_variant_ids.txt)" >> variant_list_build.log
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Process 5: Compute pairwise LD (r²) between all region variants and the
//   target locus (params.target_loci).
//   Runs once per QC stage; output named ld_<qc_stage>.vcor so the plotting
//   script can match each panel to its own LD reference.
//   Uses plink2 alpha6 (/home/b/b37974/plink2_alpha6/plink2) which supports
//   --r2-unphased --ld-snp.  BIM IDs are chr-prefixed so target_loci is used
//   directly without any prefix conversion.
// ─────────────────────────────────────────────────────────────────────────────
process COMPUTE_LD_WITH_TARGET {

    tag "LD | ${qc_stage} | ${params.target_loci}"

    executor 'slurm'
    queue    'gr10478b'
    time     '36h'

    publishDir "${params.results_dir}/02_regional_plot", mode: 'symlink'

    input:
    tuple val(qc_stage), path(bed), path(bim), path(fam)

    output:
    path "ld_${qc_stage}.vcor", emit: ld
    path "ld_${qc_stage}.log",  emit: log

    script:
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}

    # BIM variant IDs are chr-prefixed (e.g. chr16:53887925:T:C),
    # so params.target_loci can be passed directly to --ld-snp.
    /home/b/b37974/plink2_alpha6/plink2 \\
        --bfile ${bed.baseName} \\
        --r2-unphased \\
        --ld-snp "${params.target_loci}" \\
        --ld-window-kb 99999 \\
        --ld-window-r2 0 \\
        --threads ${params.plink2_threads} \\
        --out ld_${qc_stage}
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Process 6: Region-level association analysis  (additive, one job per QC stage)
//   Runs with identical pheno/covar settings across all three QC stages so that
//   results are directly comparable.  Named by QC stage for easy side-by-side
//   inspection.
// ─────────────────────────────────────────────────────────────────────────────
process RUN_REGION_ASSOC_PLINK2 {

    tag "${qc_stage} | additive | ${params.genomic_locus}"

    executor 'slurm'
    queue    'gr10478b'
    time     '36h'

    publishDir "${params.results_dir}/01_region_assoc/${qc_stage}", mode: 'symlink'

    input:
    tuple val(qc_stage), path(bed), path(bim), path(fam), path(assoc_variant_ids)

    output:
    path "*.log",                              emit: logs
    tuple val(qc_stage), path("*.glm.*"),      emit: glm

    script:
    def pheno_file         = "${params.assoc_input_root}/${params.pheno_rel_path}"
    def covar_file         = "${params.assoc_input_root}/${params.covar_rel_path}"
    def age_na_remove_file = "${params.assoc_input_root}/${params.age_na_remove_rel_path}"
    def age_na_remove_opt  = params.include_age_z ? "--remove ${age_na_remove_file}" : ''
    def locus_tag          = params.genomic_locus.replace(':', '_')
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}

    covarTag=\$(python ${params.script_dir}/covar_tag_from_names.py --covar-names "${params.covar_names}")
    outPrefix="${params.output_prefix}.region_check.${qc_stage}.${locus_tag}.${params.covar_pc_source}.\${covarTag}"

    plink2 \\
        --bfile ${bed.baseName} \\
        --pheno ${pheno_file} \\
        --pheno-name ${params.pheno_name} \\
        --covar ${covar_file} \\
        --covar-name "${params.covar_names}" \\
        ${age_na_remove_opt} \\
        --extract ${assoc_variant_ids} \
        --glm omit-ref ${params.firth_mode} hide-covar \\
        --out \${outPrefix} \\
        --ci 0.95 \\
        --threads ${params.plink2_threads}
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Process 7: Stacked regional association plot
//   Draws one panel per QC stage (shared x-axis).
//   Points are coloured by LD r² with params.target_loci.
//   Grey points = variants not present in the LD reference.
// ─────────────────────────────────────────────────────────────────────────────
process PLOT_REGIONAL_COMPARISON {

    tag "regional_plot | ${params.target_loci}"

    executor 'slurm'
    queue    'gr10478b'
    time     '36h'

    publishDir "${params.results_dir}/02_regional_plot", mode: 'symlink'

    input:
    path all_files   // staged together: *.glm.* (all QC stages) + ld_*.vcor (per stage)

    output:
    path "regional_comparison.png", emit: plot

    script:
    def geneGtfOpt = params.gene_gtf ? "--gene-gtf ${params.gene_gtf}" : ''
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}

    python ${params.script_dir}/plot_regional_comparison.py \\
        --glm-files *.glm.* \\
        --ld-files  ld_*.vcor \\
        --target-loci "${params.target_loci}" \\
        --genomic-locus "${params.genomic_locus}" \\
        ${geneGtfOpt} \\
        --title "Regional association QC-stage comparison: ${params.genomic_locus}" \\
        --out regional_comparison
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Workflow
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    validateParams()

    def chrom        = parseGenomicLocus(params.genomic_locus).chrom
    def chrAliases   = chromAliases(chrom)
    def ref_fam_path = requiredPath("${params.post_variant_qc_plink_prefix}.fam", 'post-variant-QC reference FAM')
    def ref_fam      = channel.value(ref_fam_path)

    def pre_vcf = firstExistingPath(
        chrAliases.collect { c -> "${params.pre_genotype_qc_vcf_dir}/${c}.vqc.vcf.gz" },
        'pre-genotype-QC VCF'
    )
    def pre_vcf_tbi = firstExistingPath(
        chrAliases.collect { c -> "${params.pre_genotype_qc_vcf_dir}/${c}.vqc.vcf.gz.tbi" },
        'pre-genotype-QC VCF index'
    )

    def post_gt_bed = firstExistingPath(
        chrAliases.collect { c -> "${params.post_genotype_qc_plink_dir}/${c}.plink.bed" },
        'post-genotype-QC PLINK BED'
    )
    def post_gt_bim = firstExistingPath(
        chrAliases.collect { c -> "${params.post_genotype_qc_plink_dir}/${c}.plink.bim" },
        'post-genotype-QC PLINK BIM'
    )
    def post_gt_fam = firstExistingPath(
        chrAliases.collect { c -> "${params.post_genotype_qc_plink_dir}/${c}.plink.fam" },
        'post-genotype-QC PLINK FAM'
    )

    def post_vt_bed = requiredPath("${params.post_variant_qc_plink_prefix}.bed", 'post-variant-QC PLINK BED')
    def post_vt_bim = requiredPath("${params.post_variant_qc_plink_prefix}.bim", 'post-variant-QC PLINK BIM')
    def post_vt_fam = requiredPath("${params.post_variant_qc_plink_prefix}.fam", 'post-variant-QC PLINK FAM')

    def popgmm_file = file(params.popgmm)
    log.info "PopGMM intersection subset enabled for ALL QC stages: ${popgmm_file}"

    if (params.only_validate) {
        log.info "Validation-only mode completed: all parameters and required input files passed checks."
        return
    }

    // ── Stage 1: extract target region from each QC stage ───────────────────
    pre_qc_out     = EXTRACT_REGION_PRE_GENOTYPE_QC(
        channel.value(pre_vcf),
        channel.value(pre_vcf_tbi),
        ref_fam
    )

    post_gt_qc_out = EXTRACT_REGION_POST_GENOTYPE_QC(
        channel.value(post_gt_bed),
        channel.value(post_gt_bim),
        channel.value(post_gt_fam),
        ref_fam
    )

    post_vt_qc_out = EXTRACT_REGION_POST_VARIANT_QC(
        channel.value(post_vt_bed),
        channel.value(post_vt_bim),
        channel.value(post_vt_fam)
    )

    // ── Stage 2: tag each PLINK tuple with its QC stage label ───────────────
    pre_qc_ch     = pre_qc_out.plink.map     { bed, bim, fam -> tuple("pre_genotype_qc",  bed, bim, fam) }
    post_gt_qc_ch = post_gt_qc_out.plink.map { bed, bim, fam -> tuple("post_genotype_qc", bed, bim, fam) }
    post_vt_qc_ch = post_vt_qc_out.plink.map { bed, bim, fam -> tuple("post_variant_qc",  bed, bim, fam) }

    // ── Stage 3: mandatory sample subset by intersection(ref_fam, popgmm) ──
    base_plink_ch = pre_qc_ch
        .mix(post_gt_qc_ch)
        .mix(post_vt_qc_ch)

    popgmm_subset_input_ch = base_plink_ch.map { qc_stage, bed, bim, fam ->
        tuple(qc_stage, bed, bim, fam, ref_fam_path, popgmm_file)
    }
    analysis_plink_ch = SUBSET_REGION_BY_POPGMM(popgmm_subset_input_ch).plink

    // ── Stage 4: build one shared variant list from pre-genotype-QC ────────
    pre_assoc_source_ch = analysis_plink_ch
        .filter { qc_stage, _bed, _bim, _fam -> qc_stage == 'pre_genotype_qc' }
        .map { _qc_stage, bed, bim, fam -> tuple(bed, bim, fam) }

    variant_ids_ch = BUILD_ASSOC_VARIANT_LIST(pre_assoc_source_ch).variant_ids

    // ── Stage 5: run additive association for all QC stages ────────────────
    assoc_input_ch = analysis_plink_ch

    assoc_with_variant_ch = assoc_input_ch
        .combine(variant_ids_ch)
        .map { t -> tuple(t[0], t[1], t[2], t[3], t[4]) }

    assoc_results = RUN_REGION_ASSOC_PLINK2(assoc_with_variant_ch)

    // ── Stage 6: compute LD with target locus for each QC stage ─────────────
    ld_input_ch = analysis_plink_ch

    ld_results = COMPUTE_LD_WITH_TARGET(ld_input_ch)

    // ── Stage 7: collect GLM + LD files together, then plot ────────────────
    // GLM and LD files are mixed into one collected channel so
    // PLOT_REGIONAL_COMPARISON receives a single staged directory.
    // The Python script separates them via globs: *.glm.* vs ld_*.vcor.
    all_files_ch = assoc_results.glm
        .flatMap { t -> t[1] instanceof List ? t[1] : [t[1]] }
        .mix(ld_results.ld)
        .collect()

    PLOT_REGIONAL_COMPARISON(all_files_ch)
}


