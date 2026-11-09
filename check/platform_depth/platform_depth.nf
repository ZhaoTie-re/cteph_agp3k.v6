nextflow.enable.dsl = 2

/* ============================================================================
 * check/platform_depth — read down-sampling and platform sensitivity of every signal
 *
 * Every control was sequenced on HiSeq X and every case on another platform, and
 * depth differs by platform. Two checks of the 16 GWAS loci, the Tier 1-2 genes and
 * the HLA signals, in the three nested cohorts (docs/METHODS.md):
 *   read down-sampling  the 89 NovaSeq and DNBSeq-G400RS 30x cases reduced to the
 *                       HiSeq X controls' median effective depth, re-called and
 *                       re-typed, then re-tested (baseline f = 1 + R replicates);
 *   platform sensitivity  heterogeneity across case platforms, by-platform
 *                       meta-analysis against the pooled result, leave-one-platform-out.
 * The report stage adds, per signal, the SNP-array evidence of check/array_crosscheck
 * (read only): the same people genotyped by a second technology; the control-source and
 * HLA reference-panel checks; the follow-up counts behind the flags (13.followup); and the
 * evidence matrix of every signal x every check (14.evidence), which flags and never removes.
 *
 * --run selects the stage(s), comma-separated (inputs always runs; cached after the first):
 *   inputs        00 sample table · 02 targets · 01 depth and fractions
 *   gates         G1-G4: the re-test machinery reproduces the published results on the
 *                 ORIGINAL genotypes, matrices and null models before anything is compared
 *   pilot         the pilot arms (--run_label pilot; their cost sets R, scripts/choose_R.py)
 *   downsampling  every arm: reads -> calls -> swap, HLA re-typing, re-tests; G6-G9
 *   platform      platform and control-source subsets: null fits, re-tests, statistics; G5
 *   report        summaries, array evidence, control sources, follow-up, evidence matrix,
 *                 figures, from the published results (cache false)
 * Every process, its stage, outputs and the checks it feeds: docs/PIPELINE.md (written from
 * the registry scripts/checks.py).
 *
 * Run with `bash run.sh <stage>` (README): preflight (no other driver, the frozen pieces of
 * docs/FROZEN.manifest.tsv, nextflow lint, the queue), the next run tag, a detached launch
 * resumed on the explicit session UUID, and a guard that stops the driver when a task that
 * must be cached is submitted. Never `nextflow run --help` (it runs) or -preview.
 * ==========================================================================*/

// ── locations ────────────────────────────────────────────────────────────────
params.project_dir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6'
params.script_dir  = "${projectDir}/scripts"
params.out_dir     = "${projectDir}/results"

params.model_inputs = "${params.project_dir}/wgs.auto.par/results/12_model_inputs"
params.RandomPrefix = 'cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.random_model'
params.tuning_rv    = "${params.project_dir}/tuning.rv/results"
params.assoc_saige  = "${params.project_dir}/analysis/assoc_saige/results"
params.assoc_gene   = "${params.project_dir}/analysis/assoc_gene/results"
params.assoc_hla    = "${params.project_dir}/analysis/assoc_hla.v2"
params.hla_v2       = "${params.project_dir}/hla.typing.v2"
params.array_crosscheck = "${params.project_dir}/check/array_crosscheck/results"
params.info_xlsx    = "${params.project_dir}/info/cteph_agp3k.v6.20260507.xlsx"

// ── tools (absolute: `source activate` does not put an env first on PATH) ────
params.python        = '/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3'
params.samtools      = '/LARGE0/gr10478/b37974/software/samtools/bin/samtools'
params.plink2        = '/home/b/b37974/plink2'                    // assoc_hla.v2's build
params.plink2_test   = '/home/b/b37974/plink2_alpha6/plink2'      // assoc_gene's test-VCF build
params.rvtest_bin    = '/home/b/b37974/rvtests/executable'
params.tabix         = '/home/b/b37974/htslib-1.9/tabix'
params.saige_rscript = '/home/b/b37974/anaconda3/envs/saige/bin/Rscript'
params.saige_step1   = '/home/b/b37974/anaconda3/envs/saige/bin/step1_fitNULLGLMM.R'
params.saige_step2   = '/home/b/b37974/anaconda3/envs/saige/bin/step2_SPAtests.R'
params.bcftools      = '/home/b/b37974/bcftools/bcftools'          // the cohort pipeline's bcftools 1.21
params.bcftools_plugins = '/home/b/b37974/bcftools/plugins'
params.bgzip         = '/home/b/b37974/htslib-1.9/bgzip'
params.fasta         = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa'
// GATK 4.1.4.0: the cohort's GenomicsDBImport / GenotypeGVCFs and the newest batch's HaplotypeCaller
params.gatk_jar      = '/LARGE0/gr10478/workspace/pipeline/packages/gatk-4.1.4.0/gatk-package-4.1.4.0-local.jar'
params.java8         = '/LARGE0/gr10478/workspace/pipeline/packages/jdk8u262-b10-jre/bin/java'
params.hcHeapGb      = 6                   // HaplotypeCaller heap (2 cores = 9,142 MB)
params.jointHeapGb   = 12                  // GenomicsDBImport / GenotypeGVCFs heap (4 cores = 18,284 MB)
// GATK 4.6.2.0 image: the cohort's VariantAnnotator (AlleleFraction) and VariantFiltration
params.gatk_sif      = '/home/b/b37974/simg/gatk_latest.sif'
// HLA-HD 1.7.1 (snapshot of the install that typed the cohort) and its inputs
params.hlahd_bin     = "${projectDir}/assets/hlahd.1.7.1/bin"
params.bowtie2_dir   = '/LARGE0/gr10478/b37974/software/bowtie2-2.5.5-linux-x86_64'
params.hlahd_freq    = "${params.hla_v2}/results/00.inputs/hlahd_freq_data"
params.hlahd_split   = '/LARGE0/gr10478/b37974/software/hlahd.dictionary/HLA_gene.split.3.50.0.txt'
params.hlahd_dict    = '/LARGE0/gr10478/b37974/software/hlahd.dictionary/IMGT-HLA-3.64.0_N150'
params.hla_reads     = "${params.project_dir}/hla.typing/results/01.reads"
params.hla_validator = "${params.project_dir}/hla.typing/scripts/validate_typing.py"
params.run_label     = 'downsampling'      // the results subfolder of a re-calling run (pilot / downsampling)
params.Replicates    = 3                   // replicates per deep case (3 in the pilot; final R from its timings)

// ── model settings (as the main analyses) ────────────────────────────────────
params.Cohorts    = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
params.Strata     = ['moderate_high', 'low_moderate_high']
params.HlaClasses = ['allele', 'residue']
params.CasePlatformSlugs = ['t7', 'novaseq', 'g400rs30', 'g400rs15']   // common.py PLATFORM_SLUG, case platforms
params.ControlSources = ['nagahama', 'bbj', 'htlv1', 'acc', 'lungca']    // common.py CONTROL_SOURCES, control sources
params.GeneCovar  = 'sex,pc1_avg,pc2_avg,pc3_avg,pc4_avg,pc5_avg,pc6_avg,pc7_avg,pc8_avg,pc9_avg,pc10_avg'
params.HlaCovar   = 'SEX,PC1_AVG,PC2_AVG,PC3_AVG,PC4_AVG,PC5_AVG,PC6_AVG,PC7_AVG,PC8_AVG,PC9_AVG,PC10_AVG'

// ── stage switches and sizes ─────────────────────────────────────────────────
params.run        = 'inputs'
params.depthChunk = 40          // genomes per depth job (~1 min each)

def runs() { params.run.toString().tokenize(',')*.trim() }
def src(String name) { file("${params.script_dir}/${name}", checkIfExists: true) }
def ext(String path) { file(path, checkIfExists: true) }

// ═══════════════════════════════════════════════════════════════════════════
// INPUTS
// ═══════════════════════════════════════════════════════════════════════════

/* SAMPLE_TABLE — status, platform, cohort membership, model sample sets, CRAM + index,
 * surviving original GVCFs (scripts/sample_table.py). */
process SAMPLE_TABLE {
    time '2h'
    publishDir "${params.out_dir}/00.samples", mode: 'copy'

    input:
    path script
    path lib

    output:
    path 'sample_table.tsv', emit: table
    path 'counts.tsv'
    path 'sample_summary.json'

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --out-dir .
    """
}

/* BUILD_TARGETS — signal variants, calling intervals on the original shard grid, fixed
 * sites, null SNPs (scripts/build_targets.py). */
process BUILD_TARGETS {
    time '2h'
    publishDir "${params.out_dir}/02.targets", mode: 'copy'

    input:
    path script
    path lib

    output:
    path 'hc.intervals', emit: intervals
    path 'extract.bed', emit: extract_bed
    path 'fixed_sites.tsv', emit: fixed
    path 'fixed_positions.intervals', emit: fixed_positions
    path 'signal_variants.tsv', emit: signals
    path 'null_snps.tsv', emit: null_snps
    path 'hc_intervals.tsv'
    path 'targets_qc.json'

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --out-dir .
    """
}

/* MAKE_CHUNKS — the genomes in groups of params.depthChunk for DEPTH. */
process MAKE_CHUNKS {
    executor 'local'

    input:
    path table

    output:
    path 'chunk_*.tsv', emit: chunks

    script:
    """
    ${params.python} - <<'EOF'
import pandas as pd
d = pd.read_csv('${table}', sep='\\t', dtype=str)[['sample_id', 'cram', 'crai']]
n = ${params.depthChunk}
for i in range(0, len(d), n):
    d.iloc[i:i + n].to_csv(f'chunk_{i // n:03d}.tsv', sep='\\t', index=False)
EOF
    """
}

/* DEPTH — read length and effective depth of one chunk of genomes (scripts/depth_probe.py). */
process DEPTH {
    time '4h'
    tag "${chunk.baseName}"

    input:
    path chunk
    val intervals       // the intervals as one comma-separated value, not the file: a re-run of
                        // BUILD_TARGETS that leaves them unchanged keeps these tasks cached
    path script
    path lib

    output:
    path "depth.${chunk.baseName}.tsv", emit: depth

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    echo '${intervals}' | tr ',' '\\n' > hc.intervals
    ${params.python} ${script} --chunk ${chunk} --targets hc.intervals --samtools ${params.samtools} \\
        --out depth.${chunk.baseName}.tsv
    """
}

/* DEPTH_SUMMARY — per genome and platform, target D*, fractions f_i, local depth, the
 * read-length audit (scripts/depth_summary.py). */
process DEPTH_SUMMARY {
    time '1h'
    publishDir "${params.out_dir}/01.depth", mode: 'copy'

    input:
    path depth, stageAs: 'depth/*'
    path table
    path script
    path lib

    output:
    path 'fractions.tsv', emit: fractions
    path 'depth_per_sample.tsv', emit: per_sample
    path 'platform_depth.tsv'
    path 'local_depth.tsv'
    path 'depth_summary.json'
    path 'depth_long.tsv.gz'

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --depth depth/*.tsv --samples ${table} --out-dir .
    ( head -1 \$(ls depth/*.tsv | head -1); for f in depth/*.tsv; do tail -n +2 \$f; done ) | gzip -c > depth_long.tsv.gz
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// RE-TEST MACHINERY (used by the gates on the ORIGINAL arm and by every other arm)
// ═══════════════════════════════════════════════════════════════════════════

/* SUBSET_ORIGINAL — the published genotypes at the fixed sites, rows copied byte for byte:
 * the random-model fileset (GWAS signals + context) and tuning.rv's call set (gene variants). */
process SUBSET_ORIGINAL {
    time '2h'
    tag "${cohort}"
    publishDir { "${params.out_dir}/04.genotypes/original/${cohort}" }, mode: 'copy'

    input:
    tuple val(cohort), path(fixed)
    path bedio
    path lib

    output:
    tuple val(cohort), val('original'), path('gwas.{bed,bim,fam}'), path('gene.{bed,bim,fam}'), emit: geno

    script:
    def rnd = "${params.model_inputs}/${cohort}/genotype/random_model/${params.RandomPrefix}"
    def rv  = "${params.tuning_rv}/${cohort}/02.callset_filter/filtered"
    """
    export PYTHONDONTWRITEBYTECODE=1
    awk -F'\\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[\$i] = i; next }
                 \$c["random_${cohort}"] == 1 { print \$c["variant_id"] > "gwas.ids" }
                 \$c["rv_${cohort}"] == 1     { print \$c["variant_id"] > "gene.ids" }' ${fixed}
    ${params.python} ${bedio} subset --bfile ${rnd} --ids gwas.ids --out gwas --require-all
    ${params.python} ${bedio} subset --bfile ${rv} --ids gene.ids --out gene --require-all
    """
}

/* RETEST_GWAS — SAIGE step 2 per chromosome with the cohort's ORIGINAL null model and the
 * options of assoc_saige.nf RUN_ASSOC; imputation best_guess (as published) or mean. */
process RETEST_GWAS {
    time '3h'
    tag "${arm}:${cohort}:${imputation}"
    publishDir { "${params.out_dir}/06.retest/${arm == 'original' ? 'original' : params.run_label + '/' + arm}/${cohort}" }, mode: 'copy', pattern: 'gwas.*.tsv'

    input:
    tuple val(cohort), val(arm), path(gwas), path(gene), val(imputation)
    path bedio
    path collect
    path lib

    output:
    tuple val(cohort), val(arm), val(imputation), path("gwas.${imputation}.tsv"), path('saige/chr*.txt'), emit: res
    path 'saige/chr*.log'

    script:
    def nul = "${params.assoc_saige}/${cohort}/01.assoc/main/null"
    def imp = imputation == 'best_guess' ? '' : "--impute_method=${imputation}"
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    mkdir -p saige
    for c in \$(cut -f1 gwas.bim | sort -un); do
        awk -v c=\$c '\$1 == c { print \$2 }' gwas.bim > chr\$c.ids
        ${params.python} ${bedio} subset --bfile gwas --ids chr\$c.ids --out chr\$c --require-all
        ${params.saige_rscript} ${params.saige_step2} \\
            --bedFile=chr\$c.bed --bimFile=chr\$c.bim --famFile=chr\$c.fam \\
            --AlleleOrder=alt-first --chrom=\$c \\
            --GMMATmodelFile=${nul}.rda --varianceRatioFile=${nul}.varianceRatio.txt \\
            --SAIGEOutputFile=saige/chr\$c.txt \\
            --is_Firth_beta=FALSE --LOCO=TRUE --is_output_moreDetails=TRUE ${imp} > saige/chr\$c.log 2>&1
    done
    ${params.python} ${collect} collect --saige saige/chr*.txt --arm ${arm} --cohort ${cohort} \\
        --imputation ${imputation} --out gwas.${imputation}.tsv
    """
}

/* GATE_G1 — the ORIGINAL arm (best_guess) reproduces the published saige.assoc.tsv rows. */
process GATE_G1 {
    time '1h'
    tag "${cohort}"
    publishDir "${params.out_dir}/_gates", mode: 'copy'

    input:
    tuple val(cohort), val(arm), val(imputation), path(tsv), path(raw, stageAs: 'saige/*')
    path collect

    output:
    path "G1.${cohort}.json", emit: json

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${collect} gate --saige saige/*.txt \\
        --published ${params.assoc_saige}/${cohort}/01.assoc/main/saige.assoc.tsv \\
        --cohort ${cohort} --out G1.${cohort}.json
    """
}

/* RETEST_GENE — rvtest CMC (+cmcWald) and SKAT-O on the Tier 1-2 genes, assoc_gene's inputs and
 * options, the pre- and post-flight checks of integration/gene_signals LOVO. */
process RETEST_GENE {
    time '2h'
    tag "${arm}:${cohort}:${stratum}"
    publishDir { "${params.out_dir}/06.retest/${arm == 'original' ? 'original' : params.run_label + '/' + arm}/${cohort}" }, mode: 'copy', pattern: 'gene.*.tsv'
    publishDir "${params.out_dir}/_gates", mode: 'copy', pattern: 'G2.*.json'

    input:
    tuple val(cohort), val(arm), path(gwas), path(gene), path(orig_gene, stageAs: 'orig/*'), val(stratum)
    val genes
    path gene_retest
    path bedio
    path lib

    output:
    tuple val(cohort), val(arm), val(stratum), path("gene.${stratum}.tsv"), emit: res
    path "G2.${cohort}.${stratum}.json", optional: true, emit: gate
    path "rvtest.*.log"

    script:
    def ag   = "${params.assoc_gene}/${cohort}"
    def gate = arm == 'original' ? "--gate G2.${cohort}.${stratum}.json" : ''
    def t    = task.ext.threads ?: 1
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    export PATH=${params.rvtest_bin}:\$PATH
    ${params.python} ${gene_retest} prepare --genes ${genes} --set-dir ${ag}/02.test_inputs/${stratum} \\
        --tested ${ag}/01.map/variants.${stratum}.txt --out-set set.txt --out-ids ids.txt
    ${params.plink2_test} --bfile gene --extract ids.txt --make-bed --out tested --threads ${t}
    ${params.plink2_test} --bfile tested --export vcf id-paste=iid bgz --out tested --threads ${t}
    ${params.tabix} -f -p vcf tested.vcf.gz

    pheno=${ag}/00.prep/pheno_rvt.tsv
    covar=${ag}/00.prep/covar_rvt.tsv
    coding=${ag}/00.prep/pheno_coding.json
    head -1 \$covar | tr '\\t ' '\\n' > .covar_header
    for c in \$(echo '${params.GeneCovar}' | tr ',' ' '); do
        grep -qx "\$c" .covar_header || { echo "ERROR: covariate '\$c' is absent from \$covar" >&2; exit 1; }
    done
    pheno_name=\$(${params.python} -c "import json;print(json.load(open('\$coding'))['pheno_col_rvtest'])")
    exp_case=\$(${params.python} -c "import json;print(json.load(open('\$coding'))['n_case'])")
    exp_ctrl=\$(${params.python} -c "import json;print(json.load(open('\$coding'))['n_control'])")
    for m in cmc skato; do
        if [ "\$m" = cmc ]; then opt='--burden cmc,cmcWald'; else opt='--kernel skato'; fi
        log=rvtest.${stratum}.\$m.log
        rvtest --inVcf tested.vcf.gz --pheno \$pheno --pheno-name "\$pheno_name" \\
            --covar \$covar --covar-name ${params.GeneCovar} --setFile set.txt --out rv.\$m \\
            --impute mean --noweb --numThread ${t} \$opt > "\$log" 2>&1
        grep -aq 'cannot be found' "\$log" && { echo "ERROR: rvtest dropped a covariate" >&2; exit 1; }
        line=\$(grep -a -m1 -E 'Loaded [0-9]+ cases' "\$log" || true)
        test -n "\$line" || { echo "ERROR: rvtest reported no case/control count" >&2; exit 1; }
        got_case=\$(echo "\$line" | sed -E 's/.*Loaded ([0-9]+) cases.*/\\1/')
        got_ctrl=\$(echo "\$line" | sed -E 's/.*, ([0-9]+) controls.*/\\1/')
        got_miss=\$(echo "\$line" | sed -E 's/.*and ([0-9]+) missing.*/\\1/')
        test "\$got_case" -eq "\$exp_case" && test "\$got_ctrl" -eq "\$exp_ctrl" && test "\$got_miss" -eq 0 || {
            echo "ERROR: rvtest loaded \$line; expected \$exp_case cases, \$exp_ctrl controls" >&2; exit 1; }
        if grep -aq 'Impute missing genotype to mean (by default)' "\$log"; then
            echo "ERROR: --impute never reached rvtest's parser" >&2; exit 1
        fi
        grep -aq 'Impute missing genotype to mean' "\$log" || { echo "ERROR: no imputation reported" >&2; exit 1; }
        if [ "\$m" = skato ]; then
            grep -aq 'beta1 = 1.00, beta2 = 25.00' "\$log" || { echo "ERROR: SKAT-O weights are not Beta(1,25)" >&2; exit 1; }
        fi
    done
    ${params.python} ${gene_retest} collect --arm ${arm} --cohort ${cohort} --stratum ${stratum} --set set.txt \\
        --cmc rv.cmc.CMC.assoc --wald rv.cmc.CMCWald.assoc --skato rv.skato.SkatO.assoc \\
        --coding \$coding --pheno \$pheno --gene-scan ${ag}/04.scan/${stratum}/gene_scan.tsv \\
        --denominator ${ag}/01.map/denominator.${stratum}.tsv --bfile tested --minor-from orig/gene \\
        --out gene.${stratum}.tsv ${gate}
    """
}

/* HLA_CHAIN — hla.typing.v2 steps 01 -> 02 -> 05 -> 06, unchanged, on an arm's typing archives
 * (scripts/hla_v2_wrapper.py); the ORIGINAL arm runs on the published inputs (gate G3). */
process HLA_CHAIN {
    time '6h'
    tag "${arm}"
    publishDir { "${params.out_dir}/05.hla/${arm == 'original' ? 'original' : params.run_label + '/' + arm}" }, mode: 'copy',
               pattern: 'side/{06.assoc_inputs/**,02.evidence/call_evidence.tsv,wrapper_config.json}',
               saveAs: { f -> f.replaceFirst('^side/', '') }

    input:
    tuple val(arm), val(inputs)
    path wrapper

    output:
    tuple val(arm), path('side/05.exports/allele_dosage.tsv'), path('side/05.exports/residue_dosage.tsv'),
          path('side/06.assoc_inputs'), emit: matrices
    tuple val(arm), path('side/02.evidence/call_evidence.tsv'), emit: evidence
    path 'side/06.assoc_inputs/**'
    path 'side/wrapper_config.json'

    script:
    // (2026-10-30: re-run once after the first HLA_MATRICES wrote through a staged link into side/05.exports)
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${wrapper} run --results side --inputs ${inputs} --steps 01,02,05,06
    chmod a-w side/05.exports/*.tsv
    """
}

process GATE_G3 {
    executor 'local'
    publishDir "${params.out_dir}/_gates", mode: 'copy'

    input:
    tuple val(arm), path(allele), path(residue), path(assoc_dir, stageAs: 'side/06.assoc_inputs')
    path wrapper

    output:
    path 'G3.json'

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${wrapper} gate --results side --out G3.json
    """
}

/* RETEST_HLA — one cohort x class: assoc_hla.v2's build_markers.py -> plink2 -> SAIGE step 2 on the
 * cohort's ORIGINAL null model (mean imputation) -> collect.py with the ORIGINAL tests and thresholds. */
process RETEST_HLA {
    time '1h'
    tag "${arm}:${cohort}:${cls}"
    publishDir { "${params.out_dir}/06.retest/${arm == 'original' ? 'original' : params.run_label + '/' + arm}/${cohort}" }, mode: 'copy', pattern: '{assoc,scan_qc}.random.*.tsv'

    input:
    tuple val(arm), val(cohort), val(cls), path(matrix)

    output:
    tuple val(arm), val(cohort), val(cls), path("assoc.random.${cls}.tsv"), emit: assoc
    tuple val(arm), val(cohort), val(cls), path("${cls}.bed"), path("${cls}.bim"), path("${cls}.fam"), emit: plink
    path "scan_qc.random.${cls}.tsv"
    path "random.${cls}.log"

    script:
    def S   = "${params.assoc_hla}/scripts"
    def pre = "${params.assoc_hla}/results/all/${cohort}/00.prep"
    def mk  = "${params.assoc_hla}/results/all/${cohort}/01.markers/${cls}"
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${S}/build_markers.py --case-set all --cohort ${cohort} --class ${cls} --matrix ${matrix} \\
        --coords ${params.assoc_hla}/results/_reference/hla_gene_coords.tsv --pheno ${pre}/pheno_covar.random.tsv \\
        --pheno-col PHENO1 --fixed-keep ${pre}/fixed.keep --out-vcf markers.vcf --out-info marker_info.tsv
    ${params.plink2} --vcf markers.vcf --double-id --make-bed --threads 2 --out ${cls}
    rm -f markers.vcf
    ${params.saige_rscript} ${params.saige_step2} --bedFile=${cls}.bed --bimFile=${cls}.bim --famFile=${cls}.fam \\
        --AlleleOrder=alt-first --chrom=6 --GMMATmodelFile=${pre}/null.rda \\
        --varianceRatioFile=${pre}/null.varianceRatio.txt --SAIGEOutputFile=random.${cls}.saige.txt \\
        --is_Firth_beta=FALSE --LOCO=TRUE --is_output_moreDetails=TRUE --minMAC=0.5 --maxMissing=0.15 \\
        --impute_method=mean > random.${cls}.log 2>&1
    ${params.python} ${S}/collect.py --case-set all --cohort ${cohort} --class ${cls} --model random \\
        --engine-out random.${cls}.saige.txt --engine-log random.${cls}.log --marker-info marker_info.tsv \\
        --coding ${pre}/pheno_coding.json --n-covariates 11 --thresholds ${mk}/effective_tests.${cls}.tsv \\
        --tests ${mk}/tests.${cls}.tsv --out-assoc assoc.random.${cls}.tsv --out-qc scan_qc.random.${cls}.tsv
    """
}

/* RETEST_OMNIBUS — the residue-position omnibus (fixed model, unrelated samples) of the three
 * cohorts, assoc_hla.v2's omnibus.py. */
process RETEST_OMNIBUS {
    time '1h'
    tag "${arm}"
    publishDir { "${params.out_dir}/06.retest/${arm == 'original' ? 'original' : params.run_label + '/' + arm}/_omnibus" }, mode: 'copy', pattern: 'omnibus*.residue.tsv'

    input:
    tuple val(arm), val(cohorts), path(bed, stageAs: 'c?/*'), path(bim, stageAs: 'c?/*'), path(fam, stageAs: 'c?/*')

    output:
    tuple val(arm), path('omnibus.residue.tsv'), emit: res
    path 'omnibus_*.residue.tsv'

    script:
    def S     = "${params.assoc_hla}/scripts"
    def pheno = cohorts.collect { c -> "${params.assoc_hla}/results/all/${c}/00.prep/pheno_covar.fixed.tsv" }.join(' ')
    def keep  = cohorts.collect { c -> "${params.assoc_hla}/results/all/${c}/00.prep/fixed.keep" }.join(' ')
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${S}/omnibus.py --case-set all --class residue --cohorts ${cohorts.join(',')} \\
        --bed ${bed} --pheno-fixed ${pheno} --keep ${keep} --pheno-col PHENO1 \\
        --covar-cols ${params.HlaCovar} --plink2 ${params.plink2} --alpha 0.05 --max-signals 5 --threads 2
    """
}

process GATE_G4 {
    executor 'local'
    publishDir "${params.out_dir}/_gates", mode: 'copy'

    input:
    path assoc, stageAs: 'assoc/??/*'
    path omnibus
    path script

    output:
    path 'G4.json'

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --assoc assoc/*/assoc.random.*.tsv --omnibus ${omnibus} \\
        --published-assoc ${params.assoc_hla}/results/all/_comparison/assoc_all.tsv \\
        --published-omnibus ${params.assoc_hla}/results/all/_omnibus/omnibus.residue.tsv --out G4.json
    """
}

process GATE_SUMMARY {
    executor 'local'
    publishDir "${params.out_dir}/_gates", mode: 'copy'

    input:
    path gates, stageAs: 'g/*'

    output:
    path 'gates_G1-G4.tsv'

    script:
    """
    ${params.python} - <<'PY'
import json, glob, os
rows = ['gate\\tfile\\tpass']
for f in sorted(glob.glob('g/*.json')):
    d = json.load(open(f))
    rows.append(f"{os.path.basename(f).split('.')[0]}\\t{os.path.basename(f)}\\t{d['pass']}")
open('gates_G1-G4.tsv', 'w').write('\\n'.join(rows) + '\\n')
bad = [r for r in rows[1:] if r.endswith('False')]
if bad:
    raise SystemExit('GATES FAILED: ' + '; '.join(bad))
PY
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// DOWN-SAMPLING ARMS — reads -> re-called genotypes (discovery calling, GATK 4.1.4.0)
// ═══════════════════════════════════════════════════════════════════════════

/* T7_FRACTIONS — the fraction of every down-sampled case (scripts/t7_fractions.py): the deep cases' from
 * DEPTH_SUMMARY line for line, and one fraction for the DNBSeq-T7 cases, D* / their median (user, 2026-11-03:
 * every case platform but DNBseq-G400RS 15x at the controls' median depth). */
process T7_FRACTIONS {
    executor 'local'
    publishDir "${params.out_dir}/01.depth", mode: 'copy'

    input:
    path fractions
    path per_sample
    path script

    output:
    path 'fractions.all.tsv', emit: fractions

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --fractions ${fractions} --per-sample ${per_sample} --out fractions.all.tsv
    """
}

/* ARM_PLAN — which genome runs in which arm, fraction and seed (scripts/arm_plan.py). */
process ARM_PLAN {
    executor 'local'
    publishDir "${params.out_dir}/03.calls/${params.run_label}", mode: 'copy'

    input:
    path fractions
    path table
    val mode
    path script
    path lib

    output:
    path 'arm_plan.tsv', emit: plan
    path 'orig_gvcfs.tsv', emit: gvcfs
    path 'deep_cases.tsv', emit: deep

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --fractions ${fractions} --samples ${table} --mode ${mode} \\
        --replicates ${params.Replicates} --out-plan arm_plan.tsv --out-gvcfs orig_gvcfs.tsv --out-deep deep_cases.tsv
    """
}

/* EXTRACT_READS — the calling intervals (+1 kb) of one genome, once: a multi-region query
 * (--region-file = -M -L) so a read in two overlapping regions is written once. */
process EXTRACT_READS {
    time '2h'
    tag "${sid}"

    input:
    tuple val(sid), val(cram), val(crai)
    path bed

    output:
    tuple val(sid), path("${sid}.region.bam"), path("${sid}.region.bam.bai"), emit: bam

    script:
    def t = task.ext.threads ?: 1
    """
    set -euo pipefail
    ${params.samtools} view -@ ${t} -b --region-file ${bed} --reference ${params.fasta} \\
        -o ${sid}.region.bam '${cram}##idx##${crai}'
    ${params.samtools} index ${sid}.region.bam
    """
}

/* SUBSAMPLE — one arm of one genome: f = 1 (baseline, every read) or `--subsample f --subsample-seed r`
 * (a hash of the read name: mates kept together; the same read names in every region and in the HLA
 * reads). The sample is renamed <id>__<arm> so all arms share one joint call. Gate G8 measures the
 * achieved depth over the calling intervals against f x the depth before. */
process SUBSAMPLE {
    time '1h'
    tag "${sid}:${arm}"

    input:
    tuple val(sid), val(arm), val(fraction), val(seed), path(bam), path(bai)
    val intervals

    output:
    tuple val(sid), val(arm), path("${sid}__${arm}.bam"), path("${sid}__${arm}.bam.bai"), emit: bam
    path "${sid}__${arm}.depth.tsv", emit: depth

    script:
    def t   = task.ext.threads ?: 1
    def sub = arm == 'baseline' ? '' : "--subsample ${fraction} --subsample-seed ${seed}"
    """
    set -euo pipefail
    ${params.samtools} view -@ ${t} -b ${sub} ${bam} -o arm.bam
    ${params.samtools} view -H arm.bam | awk -F'\\t' -v OFS='\\t' -v sm='SM:${sid}__${arm}' \\
        '/^@RG/ { for (i = 2; i <= NF; i++) if (\$i ~ /^SM:/) \$i = sm } { print }' > hdr.sam
    test \$(grep -c '^@RG' hdr.sam) -ge 1
    ${params.samtools} reheader hdr.sam arm.bam > ${sid}__${arm}.bam
    rm -f arm.bam
    ${params.samtools} index ${sid}__${arm}.bam
    echo -e 'sample_id\\tarm\\tfraction\\tregion\\tdepth_before\\tdepth_after' > ${sid}__${arm}.depth.tsv
    for reg in \$(echo '${intervals}' | tr ',' ' '); do
        b=\$(${params.samtools} coverage -H -q 20 -Q 10 -r \$reg ${bam} | cut -f7)
        a=\$(${params.samtools} coverage -H -q 20 -Q 10 -r \$reg ${sid}__${arm}.bam | cut -f7)
        echo -e "${sid}\\t${arm}\\t${fraction}\\t\$reg\\t\$b\\t\$a" >> ${sid}__${arm}.depth.tsv
    done
    """
}

/* CALL_GVCF — the cohort's HaplotypeCaller command (GATK 4.1.4.0, -ERC GVCF, defaults =
 * --pcr-indel-model CONSERVATIVE, MQ >= 20, BQ >= 10) on the calling intervals. Discovery calling,
 * no forced alleles: forcing replaces reference blocks by site records and changes hom-ref GQ/DP,
 * the main source of no-calls at control depth. */
process CALL_GVCF {
    time '6h'
    tag "${sid}:${arm}"

    input:
    tuple val(sid), val(arm), path(bam), path(bai)
    path intervals

    output:
    tuple val(sid), val(arm), path("${sid}__${arm}.g.vcf.gz"), path("${sid}__${arm}.g.vcf.gz.tbi"), emit: gvcf

    script:
    """
    set -euo pipefail
    ${params.java8} -Xmx${params.hcHeapGb}g -Xms${params.hcHeapGb}g -jar ${params.gatk_jar} HaplotypeCaller \\
        -R ${params.fasta} -I ${bam} -O ${sid}__${arm}.g.vcf.gz -ERC GVCF -L ${intervals} \\
        --native-pair-hmm-threads ${task.ext.threads ?: 1}
    """
}

/* JOINT — one chromosome, every GVCF of the run in one GenomicsDB (one allele set for all arms and
 * the surviving original GVCFs; one GenomicsDB partition per calling interval, since a merged partition
 * would import the whole-chromosome original GVCFs between them): variant sites over the calling
 * intervals with -stand-call-conf 0, and
 * the fixed positions with --include-non-variant-sites (4.1.4.0 has no --force-output-intervals). */
process JOINT {
    time '12h'
    tag "${chrom}"

    input:
    tuple val(chrom), val(intervals), val(positions), path(orig, stageAs: 'gvcf/*'), val(orig_lines)
    path arms, stageAs: 'gvcf/*'
    val arm_lines

    output:
    tuple val(chrom), path("${chrom}.variants.vcf.gz"), path("${chrom}.allsites.vcf.gz"), emit: vcf

    script:
    def heap  = params.jointHeapGb
    def lines = (arm_lines + orig_lines).collect { l -> "'${l}'" }.join(' ')
    """
    set -euo pipefail
    echo '${intervals}' | tr ',' '\\n' > ${chrom}.intervals
    echo '${positions}' | tr ',' '\\n' > ${chrom}.positions.intervals
    printf '%s\\n' ${lines} | awk -F'\\t' -v OFS='\\t' '{ print \$1, "gvcf/" \$2 }' > sample_map.tsv
    ${params.java8} -Xmx${heap}g -jar ${params.gatk_jar} GenomicsDBImport --sample-name-map sample_map.tsv \\
        --genomicsdb-workspace-path gdb -L ${chrom}.intervals --batch-size 50 \\
        --reader-threads ${task.ext.threads ?: 1}
    ${params.java8} -Xmx${heap}g -jar ${params.gatk_jar} GenotypeGVCFs -R ${params.fasta} -V gendb://gdb \\
        -L ${chrom}.intervals -stand-call-conf 0 -O ${chrom}.variants.vcf.gz
    ${params.java8} -Xmx${heap}g -jar ${params.gatk_jar} GenotypeGVCFs -R ${params.fasta} -V gendb://gdb \\
        -L ${chrom}.positions.intervals --only-output-calls-starting-in-intervals --include-non-variant-sites \\
        -stand-call-conf 0 -O ${chrom}.allsites.vcf.gz
    rm -rf gdb
    """
}

/* GENOTYPE_QC — the fixed sites with the cohort's own normalisation, AlleleFraction and genotype
 * filters (scripts/genotype_qc.sh, scripts/expand_fixed.py). */
process GENOTYPE_QC {
    time '2h'
    tag "${chrom}"
    publishDir "${params.out_dir}/03.calls/${params.run_label}/per_chrom", mode: 'copy', pattern: '*.qc.vcf.gz*'

    input:
    tuple val(chrom), path(variants), path(allsites)
    path fixed
    path qc_sh
    path expand

    output:
    tuple val(chrom), path("${chrom}.qc.vcf.gz"), path("${chrom}.qc.vcf.gz.tbi"), emit: vcf

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    export PYTHON=${params.python} BCFTOOLS=${params.bcftools} BGZIP=${params.bgzip} TABIX=${params.tabix}
    export FASTA=${params.fasta} GATK_SIF=${params.gatk_sif} EXPAND=${expand}
    export BCFTOOLS_PLUGINS=${params.bcftools_plugins}
    bash ${qc_sh} ${variants} ${allsites} ${fixed} ${chrom} ${chrom}
    """
}

/* CONCAT_QC — the per-chromosome QC'd calls in one VCF. */
process CONCAT_QC {
    time '1h'
    publishDir "${params.out_dir}/03.calls/${params.run_label}", mode: 'copy'

    input:
    path vcfs, stageAs: 'qc/*'

    output:
    tuple path('arms.qc.vcf.gz'), path('arms.qc.vcf.gz.tbi'), emit: vcf

    script:
    """
    set -euo pipefail
    ls qc/*.qc.vcf.gz | sort -V > list.txt
    ${params.bcftools} concat -f list.txt -Oz -o arms.qc.vcf.gz
    ${params.tabix} -f -p vcf arms.qc.vcf.gz
    """
}

/* SWAP — one arm x cohort: the 89 deep cases' QC'd genotypes replace theirs in the ORIGINAL
 * GWAS and gene filesets (scripts/bedio.py swap); every other genome is byte-identical (gate G7). */
process SWAP {
    time '1h'
    tag "${arm}:${cohort}"
    publishDir { "${params.out_dir}/04.genotypes/${arm == 'original' ? 'original' : params.run_label + '/' + arm}/${cohort}" }, mode: 'copy', pattern: '*.json'

    input:
    tuple val(cohort), val(arm), path(gwas, stageAs: 'orig/*'), path(gene, stageAs: 'orig/*'), path(vcf), path(tbi)
    path samples
    path bedio
    path lib

    output:
    tuple val(cohort), val(arm), path('gwas.{bed,bim,fam}'), path('gene.{bed,bim,fam}'), emit: geno
    path 'swap.*.json'

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    awk -F'\\t' -v a=${arm} 'NR > 1 { print \$1 "\\t" \$1 "__" a }' ${samples} > swap.samples
    ${params.python} ${bedio} swap --bfile orig/gwas --vcf ${vcf} --samples swap.samples --out gwas --report swap.gwas.json
    ${params.python} ${bedio} swap --bfile orig/gene --vcf ${vcf} --samples swap.samples --out gene --report swap.gene.json
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// DOWN-SAMPLING ARMS — HLA reads -> HLA-HD -> hla.typing.v2 matrices
// ═══════════════════════════════════════════════════════════════════════════

/* HLA_SUBSAMPLE — the archived HLA read pairs of one genome (hla.typing 01.reads, the reads that were
 * typed), subsampled like the variant reads: same fraction, same seed, same read-name hash. */
process HLA_SUBSAMPLE {
    time '1h'
    tag "${sid}:${arm}"

    input:
    tuple val(sid), val(arm), val(fraction), val(seed)

    output:
    tuple val(sid), val(arm), path("reads/${sid}.hla.R1.fastq.gz"), path("reads/${sid}.hla.R2.fastq.gz"), emit: fq

    script:
    def r1  = "${params.hla_reads}/${sid}.hla.R1.fastq.gz"
    def r2  = "${params.hla_reads}/${sid}.hla.R2.fastq.gz"
    def sub = arm == 'baseline' ? '' : "--subsample ${fraction} --subsample-seed ${seed}"
    """
    set -euo pipefail
    mkdir -p reads
    ${params.samtools} import -1 ${r1} -2 ${r2} -O bam,level=1 -o - \\
      | ${params.samtools} view -u ${sub} - \\
      | ${params.samtools} fastq -n -1 reads/${sid}.hla.R1.fastq.gz -2 reads/${sid}.hla.R2.fastq.gz \\
            -0 /dev/null -s /dev/null -
    """
}

/* HLAHD — HLA-HD 1.7.1 from the snapshot (assets/), the typing run's command, then its validator. */
process HLAHD {
    time '8h'
    tag "${sid}:${arm}"

    input:
    tuple val(sid), val(arm), path(fq1), path(fq2)

    output:
    tuple val(sid), val(arm), path("${sid}"), emit: typed

    script:
    """
    set -euo pipefail
    export PATH=${params.hlahd_bin}:${params.bowtie2_dir}:\$PATH
    mkdir -p ${sid}
    ${params.hlahd_bin}/hlahd.sh -t ${task.ext.threads ?: 1} -m 100 -c 0.95 -f ${params.hlahd_freq} \\
        ${fq1} ${fq2} ${params.hlahd_split} ${params.hlahd_dict} ${sid} ./ > hlahd.run.log 2>&1 || true
    find ${sid} -name '*.fastq' -delete || true
    rm -rf ${sid}/mapfile
    mkdir -p ${sid}/result ${sid}/log
    mv hlahd.run.log ${sid}/hlahd.run.log
    ${params.python} ${params.hla_validator} --sample ${sid} --result-dir ${sid}/result \\
        --run-log ${sid}/hlahd.run.log --gene-split ${params.hlahd_split} --out ${sid}/typing_validation.tsv
    """
}

/* HLA_INPUTS — one arm's side 00.inputs for hla.typing.v2 (scripts/hla_v2_wrapper.py inputs). */
process HLA_INPUTS {
    time '2h'
    tag "${arm}"

    input:
    tuple val(arm), val(sids), path(typed, stageAs: 'typing/*'), path(r1, stageAs: 'reads/*'), path(r2, stageAs: 'reads/*')
    path wrapper

    output:
    tuple val(arm), path('inputs'), emit: inputs

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    printf '%s\\n' ${sids.join(' ')} > retyped.txt
    ${params.python} ${wrapper} inputs --samples retyped.txt --typing-dir typing --reads-dir reads --out inputs
    """
}

/* HLA_MATRICES — one arm x cohort x class on the ORIGINAL tested columns (scripts/hla_matrices.py). */
process HLA_MATRICES {
    time '1h'
    tag "${arm}:${cohort}:${cls}"
    publishDir { "${params.out_dir}/05.hla/${arm == 'original' ? 'original' : params.run_label + '/' + arm}/matrices/${cohort}" }, mode: 'copy'

    // The exports are staged under in/ so that the output, which has the same file name, is a new file
    // here and never a write through the staged link into the HLA_CHAIN task's own 05.exports.
    input:
    tuple val(arm), path(allele, stageAs: 'in/*'), path(residue, stageAs: 'in/*'), path(assoc_dir), val(cohort), val(cls)
    path script

    output:
    tuple val(arm), val(cohort), val(cls), path("${cls}_dosage.tsv"), emit: matrix

    script:
    def chk = arm == 'original' ? '--check' : ''
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    test ! -e ${cls}_dosage.tsv
    ${params.python} ${script} --exports in --cohort ${cohort} --class ${cls} --out ${cls}_dosage.tsv ${chk}
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// PLATFORM SENSITIVITY — strata and leave-one-platform-out, all controls shared
// ═══════════════════════════════════════════════════════════════════════════

/* MAKE_SUBSETS — one cohort: the pooled prep (assoc_hla.v2 prep_pheno_cov.py --case-set all, must equal
 * the published 00.prep), 4 platform strata and 4 leave-one-platform-out sets (scripts/make_subsets.py). */
process MAKE_SUBSETS {
    time '1h'
    tag "${cohort}"
    publishDir { "${params.out_dir}/08.subsets/${cohort}" }, mode: 'copy', saveAs: { f -> f.replaceFirst('^subsets/', '') }

    input:
    tuple val(cohort), path(table)
    path script
    path lib

    output:
    tuple val(cohort), path('subsets'), emit: dir
    path 'subsets/*'

    script:
    def S   = "${params.assoc_hla}/scripts"
    def mi  = "${params.model_inputs}/${cohort}"
    def fam = "${mi}/genotype/fixed_model/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.maf_ge_threshold.fam"
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${S}/prep_pheno_cov.py --cohort ${cohort} --case-set all \\
        --case-sets ${params.assoc_hla}/results/_reference/case_sets.tsv --pheno ${mi}/phenotype/pheno.tsv \\
        --pheno-col PHENO1 --covar ${mi}/covariates/bbj_mainland_pc.sex.tsv --covar-cols ${params.HlaCovar} \\
        --fixed-fam ${fam} --out-fixed pooled.fixed.tsv --out-random pooled.random.tsv --out-keep pooled.keep \\
        --out-json pooled.json
    ${params.python} ${script} --cohort ${cohort} --samples ${table} --random pooled.random.tsv --fixed pooled.fixed.tsv \\
        --keep pooled.keep --published-prep ${params.assoc_hla}/results/all/${cohort}/00.prep --out-dir subsets
    """
}

/* MAKE_CTRL_SUBSETS — one cohort: every case with one source's controls (ctrlonly_<source>) or with all
 * sources but one (ctrlwithout_<source>); the rows of the pooled platform subset, unchanged
 * (scripts/control_subsets.py; sources = the info sheet's Project column). The folder keeps the name
 * 'subsets' so PLATFORM_OMNIBUS reads it as it reads the platform subsets. */
process MAKE_CTRL_SUBSETS {
    time '1h'
    tag "${cohort}"
    publishDir { "${params.out_dir}/08.subsets/${cohort}/control_sources" }, mode: 'copy'

    input:
    tuple val(cohort), path(platform, stageAs: 'platform')
    path info
    path script
    path lib

    output:
    tuple val(cohort), path('subsets'), emit: dir

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --cohort ${cohort} --subsets platform --info ${info} --out-dir subsets
    """
}

/* FIT_NULL — SAIGE step 1 on one subset, the genome-wide scan's settings and GRM markers (a copy of
 * assoc_hla.v2 FIT_NULL: its grm.prune.in on the random-model genotypes, the subset's samples). */
process FIT_NULL {
    time '12h'
    tag "${cohort}:${subset}"
    publishDir { "${params.out_dir}/08.subsets/${cohort}/null/${subset}" }, mode: 'copy', pattern: 'null.{varianceRatio.txt,fit.log}'

    input:
    tuple val(cohort), val(subset), path(pheno_random)

    output:
    tuple val(cohort), val(subset), path('null.rda'), path('null.varianceRatio.txt'), emit: null_model
    path 'null.fit.log'

    script:
    def rnd = "${params.model_inputs}/${cohort}/genotype/random_model/${params.RandomPrefix}"
    def t   = task.ext.threads ?: 1
    """
    set -euo pipefail
    tail -n +2 ${pheno_random} | cut -f1,2 > samples.keep
    ${params.plink2} --bed ${rnd}.bed --bim ${rnd}.bim --fam ${rnd}.fam \\
        --extract ${params.assoc_saige}/${cohort}/00.prep/grm.prune.in --keep samples.keep --make-bed \\
        --out grm_markers --threads ${t}
    sed '1s/^#FID/FID/' ${pheno_random} > pheno_covar.null.tsv
    ${params.saige_rscript} ${params.saige_step1} \\
        --plinkFile=grm_markers --phenoFile=pheno_covar.null.tsv --phenoCol=PHENO1 \\
        --covarColList=${params.HlaCovar} --sampleIDColinphenoFile=IID --traitType=binary --outputPrefix=null \\
        --nThreads=${t} --isDiagofKinSetAsOne=True --numRandomMarkerforVarianceRatio=200 \\
        --skipVarianceRatioEstimation=FALSE --useSparseGRMtoFitNULL=FALSE --LOCO=TRUE \\
        --IsOverwriteVarianceRatioFile=TRUE --isCovariateOffset=FALSE 2>&1 | tee null.fit.log
    rm -f grm_markers.bed grm_markers.bim grm_markers.fam
    """
}

/* PLATFORM_GENO — one cohort: the GWAS signal fileset joined with the null SNPs (random-model rows,
 * byte copies), for every subset's step 2. */
process PLATFORM_GENO {
    time '2h'
    tag "${cohort}"

    input:
    tuple val(cohort), path(fixed), path(nulls)
    path bedio
    path lib

    output:
    tuple val(cohort), path('pmark.{bed,bim,fam}'), emit: geno

    script:
    def rnd = "${params.model_inputs}/${cohort}/genotype/random_model/${params.RandomPrefix}"
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    awk -F'\\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[\$i] = i; next }
                 \$c["random_${cohort}"] == 1 { print \$c["variant_id"] }' ${fixed} > ids.txt
    awk -F'\\t' 'NR > 1 { print \$1 }' ${nulls} >> ids.txt
    sort -u ids.txt > ids.u.txt
    ${params.python} ${bedio} subset --bfile ${rnd} --ids ids.u.txt --out pmark
    """
}

/* STEP2_SUBSET — SAIGE step 2 of one subset's null model on the GWAS signals + null SNPs (per
 * chromosome; best_guess as the published scan) and the HLA alleles and residues (mean, as assoc_hla.v2). */
process STEP2_SUBSET {
    time '4h'
    tag "${cohort}:${subset}"
    publishDir { "${params.out_dir}/09.platform_assoc/${cohort}" }, mode: 'copy', pattern: '*.saige.tsv'

    input:
    tuple val(cohort), val(subset), path(rda), path(vr), path(pmark), path(allele, stageAs: 'hla/*'),
          path(residue, stageAs: 'hla/*')
    path bedio
    path collect
    path lib

    output:
    tuple val(cohort), val(subset), path("${subset}.saige.tsv"), emit: res

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    mkdir -p out
    for c in \$(cut -f1 pmark.bim | sort -un); do
        awk -v c=\$c '\$1 == c { print \$2 }' pmark.bim > chr\$c.ids
        ${params.python} ${bedio} subset --bfile pmark --ids chr\$c.ids --out chr\$c --require-all
        ${params.saige_rscript} ${params.saige_step2} --bedFile=chr\$c.bed --bimFile=chr\$c.bim --famFile=chr\$c.fam \\
            --AlleleOrder=alt-first --chrom=\$c --GMMATmodelFile=${rda} --varianceRatioFile=${vr} \\
            --SAIGEOutputFile=out/chr\$c.txt --is_Firth_beta=FALSE --LOCO=TRUE --is_output_moreDetails=TRUE \\
            > out/chr\$c.log 2>&1
        rm -f chr\$c.bed chr\$c.bim chr\$c.fam
    done
    for k in allele residue; do
        ${params.saige_rscript} ${params.saige_step2} --bedFile=hla/\$k.bed --bimFile=hla/\$k.bim --famFile=hla/\$k.fam \\
            --AlleleOrder=alt-first --chrom=6 --GMMATmodelFile=${rda} --varianceRatioFile=${vr} \\
            --SAIGEOutputFile=out/hla_\$k.txt --is_Firth_beta=FALSE --LOCO=TRUE --is_output_moreDetails=TRUE \\
            --minMAC=0.5 --maxMissing=0.15 --impute_method=mean > out/hla_\$k.log 2>&1
    done
    ${params.python} ${collect} collect --saige out/*.txt --arm original --cohort ${cohort} --subset ${subset} \\
        --imputation mixed --out ${subset}.saige.tsv
    """
}

/* PLATFORM_GENE — rvtest CMC / SKAT-O of the Tier 1-2 genes on one subset (--peopleIncludeFile). */
process PLATFORM_GENE {
    time '2h'
    tag "${cohort}:${subset}:${stratum}"
    publishDir { "${params.out_dir}/09.platform_assoc/${cohort}" }, mode: 'copy', pattern: 'gene.*.tsv'

    input:
    tuple val(cohort), val(subset), path(incl), path(gene), val(stratum)
    val genes
    path gene_retest
    path bedio
    path lib

    output:
    tuple val(cohort), val(subset), val(stratum), path("gene.${subset}.${stratum}.tsv"), emit: res

    script:
    def ag = "${params.assoc_gene}/${cohort}"
    def t  = task.ext.threads ?: 1
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    export PATH=${params.rvtest_bin}:\$PATH
    ${params.python} ${gene_retest} prepare --genes ${genes} --set-dir ${ag}/02.test_inputs/${stratum} \\
        --tested ${ag}/01.map/variants.${stratum}.txt --out-set set.txt --out-ids ids.txt
    ${params.plink2_test} --bfile gene --extract ids.txt --make-bed --out tested --threads ${t}
    ${params.plink2_test} --bfile tested --export vcf id-paste=iid bgz --out tested --threads ${t}
    ${params.tabix} -f -p vcf tested.vcf.gz
    pheno=${ag}/00.prep/pheno_rvt.tsv
    covar=${ag}/00.prep/covar_rvt.tsv
    coding=${ag}/00.prep/pheno_coding.json
    pheno_name=\$(${params.python} -c "import json;print(json.load(open('\$coding'))['pheno_col_rvtest'])")
    for m in cmc skato; do
        if [ "\$m" = cmc ]; then opt='--burden cmc,cmcWald'; else opt='--kernel skato'; fi
        log=rvtest.\$m.log
        rvtest --inVcf tested.vcf.gz --pheno \$pheno --pheno-name "\$pheno_name" --peopleIncludeFile ${incl} \\
            --covar \$covar --covar-name ${params.GeneCovar} --setFile set.txt --out rv.\$m \\
            --impute mean --noweb --numThread ${t} \$opt > "\$log" 2>&1
        grep -aq 'cannot be found' "\$log" && { echo "ERROR: rvtest dropped a covariate" >&2; exit 1; }
        grep -aq 'Impute missing genotype to mean (by default)' "\$log" && { echo "ERROR: --impute not parsed" >&2; exit 1; }
        grep -aq 'Impute missing genotype to mean' "\$log" || { echo "ERROR: no imputation reported" >&2; exit 1; }
        if [ "\$m" = skato ]; then
            grep -aq 'beta1 = 1.00, beta2 = 25.00' "\$log" || { echo "ERROR: SKAT-O weights are not Beta(1,25)" >&2; exit 1; }
        fi
    done
    ${params.python} ${gene_retest} collect --arm original --subset ${subset} --cohort ${cohort} --stratum ${stratum} \\
        --set set.txt --cmc rv.cmc.CMC.assoc --wald rv.cmc.CMCWald.assoc --skato rv.skato.SkatO.assoc \\
        --coding \$coding --pheno \$pheno --gene-scan ${ag}/04.scan/${stratum}/gene_scan.tsv \\
        --denominator ${ag}/01.map/denominator.${stratum}.tsv --bfile tested --minor-from gene \\
        --include ${incl} --out gene.${subset}.${stratum}.tsv
    """
}

/* CARRIER_FIRTH — one cohort: the carrier OR of each gene x variant set in every subset (descriptive). */
process CARRIER_FIRTH {
    time '1h'
    tag "${cohort}"
    publishDir { "${params.out_dir}/09.platform_assoc/${cohort}" }, mode: 'copy', pattern: 'carriers_firth.tsv'

    input:
    tuple val(cohort), path(gene), path(subsets), val(genes)
    path gene_retest
    path firth
    path bedio
    path lib

    output:
    tuple val(cohort), path('carriers_firth.tsv'), path('carriers.{bed,bim,fam}'), emit: res

    script:
    def ag   = "${params.assoc_gene}/${cohort}"
    def subs = (['pooled'] + params.CasePlatformSlugs.collect { s -> "only_${s}" } +
                params.CasePlatformSlugs.collect { s -> "without_${s}" })
    def inc  = subs.collect { s -> "${s}=${subsets}/${s}.rv_include.txt" }.join(' ')
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    for s in ${params.Strata.join(' ')}; do
        ${params.python} ${gene_retest} prepare --genes ${genes} --set-dir ${ag}/02.test_inputs/\$s \\
            --tested ${ag}/01.map/variants.\$s.txt --out-set set.\$s.txt --out-ids ids.\$s.txt
    done
    ${params.python} ${firth} --bfile gene --sets ${params.Strata.collect { s -> "${s}=set.${s}.txt" }.join(' ')} \\
        --pheno ${subsets}/pooled.fixed.tsv --include ${inc} --plink2 ${params.plink2} \\
        --covar-cols SEX,PC1_AVG-PC10_AVG --cohort ${cohort} --out carriers_firth.tsv
    """
}

/* PLATFORM_OMNIBUS — the residue-position omnibus of the three cohorts on one subset. */
process PLATFORM_OMNIBUS {
    time '1h'
    tag "${subset}"
    publishDir { "${params.out_dir}/09.platform_assoc/_omnibus" }, mode: 'copy', pattern: 'omnibus.*.tsv'

    input:
    tuple val(subset), val(cohorts), path(bed, stageAs: 'c?/*'), path(bim, stageAs: 'c?/*'), path(fam, stageAs: 'c?/*'),
          path(dirs, stageAs: 'sub?/*')

    output:
    tuple val(subset), path("omnibus.${subset}.tsv"), emit: res

    script:
    def S     = "${params.assoc_hla}/scripts"
    def pheno = (1..cohorts.size()).collect { i -> "sub${i}/subsets/${subset}.fixed.tsv" }.join(' ')
    def keep  = (1..cohorts.size()).collect { i -> "sub${i}/subsets/${subset}.keep" }.join(' ')
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${S}/omnibus.py --case-set all --class residue --cohorts ${cohorts.join(',')} \\
        --bed ${bed} --pheno-fixed ${pheno} --keep ${keep} --pheno-col PHENO1 \\
        --covar-cols ${params.HlaCovar} --plink2 ${params.plink2} --alpha 0.05 --max-signals 1 --threads 2
    mv omnibus.residue.tsv omnibus.${subset}.tsv
    """
}

/* PLATFORM_STATS — one cohort: heterogeneity (Q, I2, empirical P), FE / RE meta vs pooled, LOPO with
 * expected z and difference tests, the shared-control correlation, null calibration, case-only checks. */
process PLATFORM_STATS {
    time '4h'
    tag "${cohort}"
    publishDir { "${params.out_dir}/10.platform/${cohort}" }, mode: 'copy'

    input:
    tuple val(cohort), val(subsets), path(saige, stageAs: 'saige/*'), path(subdir), path(pmark),
          path(carriers), path(carriers_geno)
    path nulls
    path signals
    path table
    path script
    path bedio
    path lib

    output:
    tuple val(cohort), path('platform.*'), emit: res

    script:
    // a staged input already carries its stageAs folder (saige/<name>)
    def specs = [subsets, saige].transpose().collect { s, f -> "${s}=${f}" }.join(' ')
    def mats  = params.HlaClasses.collect { k -> "${params.hla_v2}/results/06.assoc_inputs/${cohort}/${k}_dosage.tsv" }.join(' ')
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    ${params.python} ${script} --cohort ${cohort} --counts ${subdir}/subset_counts.tsv --saige ${specs} \\
        --null ${nulls} --signals ${signals} --samples ${table} --pheno-fixed ${subdir}/pooled.fixed.tsv \\
        --gwas-bfile pmark --hla-matrices ${mats} --carriers-bfile carriers --covar-cols ${params.HlaCovar} \\
        --out-prefix platform
    """
}

/* GATE_G5 — the pooled re-fit reproduces the published null model: variance ratio and step-2 z of
 * the GWAS signals within 1 % (the variance-ratio markers are drawn at random by SAIGE). */
process GATE_G5 {
    executor 'local'
    tag "${cohort}"
    publishDir "${params.out_dir}/_gates", mode: 'copy'

    input:
    tuple val(cohort), path(vr), path(refit), path(original)

    output:
    path "G5.${cohort}.json"

    script:
    """
    ${params.python} - <<'PY'
import json
import numpy as np
import pandas as pd
pub = float(open('${params.assoc_saige}/${cohort}/01.assoc/main/null.varianceRatio.txt').read().split()[0])
new = float(open('${vr}').read().split()[0])
a = pd.read_csv('${refit}', sep='\\t', dtype={'variant_id': str}).set_index('variant_id')
b = pd.read_csv('${original}', sep='\\t', dtype={'variant_id': str}).set_index('variant_id')
ids = a.index.intersection(b.index)
rel = float(np.max(np.abs(a.loc[ids, 'z'] - b.loc[ids, 'z']) / np.maximum(np.abs(b.loc[ids, 'z']), 1e-12)))
med = float(np.median(np.abs(a.loc[ids, 'z'] - b.loc[ids, 'z']) / np.maximum(np.abs(b.loc[ids, 'z']), 1e-12)))
rep = {'cohort': '${cohort}', 'variance_ratio_published': pub, 'variance_ratio_refit': new,
       'variance_ratio_rel_diff': abs(new - pub) / pub, 'variants': int(len(ids)),
       'z_max_rel_diff': rel, 'z_median_rel_diff': med}
rep['pass'] = bool(rep['variance_ratio_rel_diff'] < 0.01 and med < 0.01)
json.dump(rep, open('G5.${cohort}.json', 'w'), indent=2)
print(rep)
PY
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// CONCORDANCE (user request 2026-10-30): genotypes, QC-rule failures, HLA calls
// ═══════════════════════════════════════════════════════════════════════════

/* CONCORDANCE — variant genotypes at the fixed sites: re-processing, joint genotyping, the
 * HaplotypeCaller re-run, the depth effect, the total change (scripts/concordance.py variants). */
process CONCORDANCE {
    time '2h'
    publishDir "${params.out_dir}/07.downsampling/${params.run_label}/concordance", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    tuple path(gwas, stageAs: 'orig_gwas/*'), path(gene, stageAs: 'orig_gene/*')
    path fixed
    path plan
    path script
    path bedio

    output:
    path 'concordance.summary.tsv', emit: summary
    path 'concordance.{long.tsv.gz,samples.tsv,replicates.tsv}'

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} variants --vcf ${vcf} --original orig_gwas/gwas orig_gene/gene --fixed ${fixed} \\
        --plan ${plan} --out-long concordance.long.tsv --out-summary concordance.summary.tsv \\
        --out-samples concordance.samples.tsv --out-replicates concordance.replicates.tsv
    gzip -f concordance.long.tsv
    """
}

/* QC_REASONS — the share of genotypes each of the cohort's QC rules no-calls: the arms (FT written
 * by VariantFiltration) and the cohort's own pre-QC genotypes by group (scripts/concordance.py reasons). */
process QC_REASONS {
    time '3h'
    publishDir "${params.out_dir}/07.downsampling/${params.run_label}/qc_reasons", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    path fixed
    path plan
    path table
    path script
    path bedio

    output:
    path 'qc_reasons.*.tsv'

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} reasons --vcf ${vcf} --fixed ${fixed} --plan ${plan} --samples ${table} \\
        --cohort-vcf ${params.project_dir}/wgs.auto.par/results/03_annotate_af_norm_gt/@CHROM@.vqc.af.gtnorm.vcf.gz \\
        --bcftools ${params.bcftools} --out-samples qc_reasons.samples.tsv --out-summary qc_reasons.summary.tsv
    """
}

/* HLA_CONCORDANCE — the deep cases' HLA calls and tested dosages, arm against archive and baseline
 * (scripts/hla_concordance.py). */
process HLA_CONCORDANCE {
    time '1h'
    publishDir "${params.out_dir}/07.downsampling/${params.run_label}/hla", mode: 'copy'

    input:
    val arm_specs       // <arm>=<absolute path of the arm's call_evidence.tsv>
    val mat_specs       // <arm>:<cohort>:<class>=<absolute path of the arm's matrix>
    path deep
    path table
    path script

    output:
    path 'hla.calls_long.tsv', emit: calls
    path 'hla.{summary,read_support,dosage}.tsv'

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --archive ${params.hla_v2}/results/02.evidence/call_evidence.tsv \\
        --arms ${arm_specs.join(' ')} --deep ${deep} --samples ${table} --matrices ${mat_specs.join(' ')} \\
        --original-matrices ${params.hla_v2}/results/06.assoc_inputs --out-prefix hla
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// REPORT — reads the published results (run after the stages above); never cached
// ═══════════════════════════════════════════════════════════════════════════

/* DS_SUMMARY — every signal across arms, and the exposure of each signal in the deep cases. */
process DS_SUMMARY {
    cache false
    time '1h'
    publishDir "${params.out_dir}/07.downsampling/${params.run_label}/summary", mode: 'copy'

    input:
    path script
    path libs       // the modules the report scripts import (from scripts/ itself; staged to declare them)

    output:
    path 'ds.*', emit: tables

    script:
    def r = params.out_dir
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --retest ${r}/06.retest --label ${params.run_label} \\
        --signals ${r}/02.targets/signal_variants.tsv --deep ${r}/03.calls/${params.run_label}/deep_cases.tsv \\
        --samples ${r}/00.samples/sample_table.tsv --orig-geno ${r}/04.genotypes/original --out-prefix ds
    """
}

/* PLATFORM_REPORT — the platform tables across cohorts (scripts/platform_report.py). */
process PLATFORM_REPORT {
    cache false
    time '1h'
    publishDir "${params.out_dir}/10.platform/_summary", mode: 'copy'

    input:
    path script
    path libs

    output:
    path 'platform.*.tsv', emit: tables

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --results ${params.out_dir} --out-dir .
    """
}

/* ARRAY_COHORTS, ARRAY_EVIDENCE — WGS against SNP arrays in the same people, in each cohort. ARRAY_COHORTS (cached):
 * scripts/array_cohorts.py re-runs check/array_crosscheck's steps on its inputs (read only) for each cohort's
 * unrelated DNBSeq-T7 cases and HiSeq X controls with arrays (gate G10: full mainland reproduces array_crosscheck's
 * tables), then scripts/array_hibag.py imputes DQB1 with HIBAG. ARRAY_EVIDENCE (rerun, seconds):
 * scripts/array_evidence.py reads one marker per signal (typed, imputed, LD proxy; DQB1:57:V with HIBAG) and the
 * gene-driver carriers on the imputed arrays (full mainland). Work files keyed by sample stay in the task folder. */
process ARRAY_COHORTS {
    // cached (user 2026-11-05: the report reruns no longer repeat it): its scripts are staged, and what it reads
    // (check/array_crosscheck, the QC'd and imputed arrays, the model inputs) is frozen (docs/FROZEN.manifest.tsv)
    time '3h'
    clusterOptions '--rsc p=1:t=4:c=4:m=18284M'
    publishDir "${params.out_dir}/11.array", mode: 'copy'

    input:
    path libs       // common.py, array_cohorts.py, array_hibag.py and array_hibag.R

    output:
    path 'cohorts', emit: cohorts

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} array_cohorts.py --array ${params.array_crosscheck} --plink2 ${params.plink2} \\
        --plink1 /home/b/b37974/plink --threads 4 --memory 14000 --work array_work --out cohorts
    # DQB1:57:V on the arrays: HIBAG with its published Asian model (assets/hibag; user 2026-11-05)
    ${params.python} array_hibag.py --array ${params.array_crosscheck} --plink2 ${params.plink2} \\
        --cohorts cohorts --work hibag_work --threads 4
    """
}

process ARRAY_EVIDENCE {
    cache false
    time '30m'
    publishDir "${params.out_dir}/11.array", mode: 'copy'

    input:
    path script
    path libs
    path cohorts    // ARRAY_COHORTS

    output:
    path 'array.*', emit: tables

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --cohorts ${cohorts} \\
        --signals ${params.out_dir}/02.targets/signal_variants.tsv --out-prefix array
    """
}

/* CONTROL_REPORT — control-source sensitivity of every signal (scripts/control_report.py): each source's
 * controls alone and each left out, from the published per-subset SAIGE, rvtest and omnibus results; the
 * source estimates share the cases (Lin & Sullivan's correlation with the roles swapped). */
process CONTROL_REPORT {
    cache false
    time '1h'
    publishDir "${params.out_dir}/12.controls", mode: 'copy'

    input:
    path script
    path libs       // with platform_stats.py

    output:
    path 'controls.*', emit: tables

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --results ${params.out_dir} --out-dir .
    """
}

/* HLA_REF_SHIFT — are the cases shifted toward ToMMo at HLA (scripts/hla_reference_shift.py)? The slope of
 * the case-control differences on the ToMMo-control differences over class I and DP: frequencies,
 * published (PC-adjusted) log ORs, and log ORs given DQB1:57:V; DQB1:57:V predicted from each. */
process HLA_REF_SHIFT {
    cache false
    time '2h'
    publishDir "${params.out_dir}/12.controls", mode: 'copy'

    input:
    path script
    path libs

    output:
    path 'hla_reference_shift.*', emit: tables

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --exports ${params.hla_v2}/results/05.exports \\
        --conditional ${params.assoc_hla}/results/all --results ${params.out_dir} --out-prefix hla_reference_shift
    """
}

/* FOLLOWUP — what lies behind the flags, in each cohort (scripts/signal_followup.py): GWAS lead frequencies by
 * case platform, registration study and control source beside ToMMo and 1KG JPT; gene carriers by group, with
 * the cohort's genotype filters re-applied to the pre-QC records (carriers set missing); the deep cases whose
 * DQB1:57:V dosage changes between arms. Counts, no new test. */
process FOLLOWUP {
    cache false
    time '3h'
    publishDir "${params.out_dir}/13.followup", mode: 'copy'

    input:
    path script
    path libs

    output:
    path 'followup.*', emit: tables

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --results ${params.out_dir} --label ${params.run_label} \\
        --cohort-vcf '${params.project_dir}/wgs.auto.par/results/03_annotate_af_norm_gt/@CHROM@.vqc.af.gtnorm.vcf.gz' \\
        --bcftools ${params.bcftools} --out-dir .
    """
}

/* EVIDENCE — every signal x every check x cohort, read by the rules of the registry (scripts/checks.py,
 * scripts/evidence_matrix.py); flags, never removes. */
process EVIDENCE {
    cache false
    time '1h'
    publishDir "${params.out_dir}/14.evidence", mode: 'copy'

    input:
    path script
    path libs
    path tables     // ds.*, array.*, controls.*, followup.* of this report run

    output:
    path 'evidence.*', emit: tables

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --results ${params.out_dir} --label ${params.run_label} --ds-summary . --array-tables . \\
        --control-tables . --followup-tables . --out-dir .
    """
}

/* FIGURES + CATALOGUE — the six figure families (scripts/figures.py) and one README each (scripts/catalogue.py). */
process FIGURES {
    cache false
    time '1h'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path tables
    path script
    path catalogue
    path libs

    output:
    path 'figures/**'

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --results ${params.out_dir} --label ${params.run_label} --out figures --ds-summary . --array-tables . --control-tables . \\
        --followup-tables . --evidence-tables .
    ${params.python} ${catalogue} --figures figures
    """
}

/* GATES_DS — G6 re-processing fidelity, G8 achieved depth, G9 HLA re-typing at f = 1 (scripts/gates_ds.py). */
process GATES_DS {
    executor 'local'
    publishDir "${params.out_dir}/_gates/${params.run_label}", mode: 'copy'

    input:
    path summary
    path depth, stageAs: 'depth/*'
    path calls
    path script

    output:
    path 'G{6,8,9}.json'

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --concordance ${summary} --depth depth/*.tsv --hla ${calls} --out-dir .
    """
}

// ═══════════════════════════════════════════════════════════════════════════
workflow {
    def lib     = src('common.py')
    def bedio   = src('bedio.py')
    def wrapper = src('hla_v2_wrapper.py')

    // ── inputs (always; cached after the first run) ──────────────────────────
    table   = SAMPLE_TABLE(src('sample_table.py'), lib).table
    targets = BUILD_TARGETS(src('build_targets.py'), lib)
    chunks  = MAKE_CHUNKS(table).chunks.flatten()
    ivals   = targets.intervals.map { f -> f.readLines()*.trim().findAll { l -> l }.join(',') }
    depth   = DEPTH(chunks, ivals, src('depth_probe.py'), lib).depth
    dsum    = DEPTH_SUMMARY(depth.collect(), table, src('depth_summary.py'), lib)

    // the Tier 1-2 genes, from the signal table (one comma-separated value)
    genes = targets.signals.map { f ->
        f.readLines().drop(1).collect { l -> l.split('\t') }.findAll { c -> c[0] == 'gene' }
            .collect { c -> c[1] }.unique().sort().join(',')
    }

    def doDs = ('pilot' in runs()) || ('downsampling' in runs())
    def doPf = 'platform' in runs()
    if (('gates' in runs()) || doDs || doPf) {

        // ── the ORIGINAL arm: published genotypes, typing and null models ────────
        orig = SUBSET_ORIGINAL(channel.fromList(params.Cohorts).combine(targets.fixed), bedio, lib).geno
        geno = orig
        chain_in = channel.of(tuple('original', "${params.hla_v2}/results/00.inputs"))
        plan  = null
        calls = null
        sub_depth = null

        // ── down-sampling arms ───────────────────────────────────────────────────
        if (doDs) {
            def mode = 'pilot' in runs() ? 'pilot' : 'full'
            fr_all = T7_FRACTIONS(dsum.fractions, dsum.per_sample, src('t7_fractions.py')).fractions
            plan = ARM_PLAN(fr_all, table, mode, src('arm_plan.py'), lib)
            rows = plan.plan.splitCsv(header: true, sep: '\t')

            reg = EXTRACT_READS(rows.map { r -> tuple(r.sample_id, r.cram, r.crai) }.unique(), targets.extract_bed).bam
            sub = SUBSAMPLE(rows.map { r -> tuple(r.sample_id, r.arm, r.fraction, r.seed) }.combine(reg, by: 0), ivals)
            sub_depth = sub.depth
            gv  = CALL_GVCF(sub.bam, targets.intervals).gvcf

            // one joint call per chromosome: every arm + the surviving original GVCFs
            chrom_iv = targets.intervals.flatMap { f ->
                f.readLines().findAll { l -> l.trim() }.groupBy { l -> l.split(':')[0] }
                    .collect { c, ls -> tuple(c, ls.join(',')) }
            }
            chrom_pos = targets.fixed_positions.flatMap { f ->
                f.readLines().findAll { l -> l.trim() }.groupBy { l -> l.split(':')[0] }
                    .collect { c, ls -> tuple(c, ls.join(',')) }
            }
            orig_g = chrom_iv.map { c, _iv -> c }
                .combine(plan.gvcfs.splitCsv(header: true, sep: '\t'))
                .map { c, r ->
                    def g = file("${r.orig_gvcf_dir}/${r.sample_id}.${c}.g.vcf.gz", checkIfExists: true)
                    tuple(c, [g, file("${g}.tbi", checkIfExists: true)], "${r.sample_id}__origgvcf\t${g.name}")
                }
                .groupTuple()
                .map { c, fs, ls -> tuple(c, fs.flatten(), ls) }
            jin = chrom_iv.join(chrom_pos).join(orig_g, remainder: true)
                .map { c, iv, ps, fs, ls -> tuple(c, iv, ps, fs ?: [], ls ?: []) }
            arm_files = gv.map { _s, _a, g, t -> [g, t] }.flatten().collect()
            arm_lines = gv.map { s, a, g, _t -> "${s}__${a}\t${g.name}" }.collect()
            jv = JOINT(jin, arm_files, arm_lines).vcf
            qc = GENOTYPE_QC(jv, targets.fixed, src('genotype_qc.sh'), src('expand_fixed.py')).vcf
            calls = CONCAT_QC(qc.map { _c, v, t -> [v, t] }.flatten().collect()).vcf

            // the deep cases' arms swapped into the ORIGINAL filesets
            deep_arms = rows.filter { r -> r.role == 'deep' }.map { r -> r.arm }.unique()
            swap_in = orig.combine(deep_arms).map { c, _o, gw, gn, a -> tuple(c, a, gw, gn) }.combine(calls)
            sw = SWAP(swap_in, plan.deep, bedio, lib).geno
            geno = orig.mix(sw)

            // HLA: the deep cases' archived HLA reads, the same fraction and seed
            hfq = HLA_SUBSAMPLE(rows.filter { r -> r.role == 'deep' }.map { r -> tuple(r.sample_id, r.arm, r.fraction, r.seed) }).fq
            typed = HLAHD(hfq).typed
            by_arm = typed.join(hfq, by: [0, 1])
                .map { s, a, t, r1, r2 -> tuple(a, s, t, r1, r2) }
                .groupTuple()
            hin = HLA_INPUTS(by_arm, wrapper).inputs
            chain_in = chain_in.mix(hin.map { a, p -> tuple(a, p.toString()) })
        }

        // ── HLA matrices: hla.typing.v2 steps 01-06 on each arm (gate G3 on the ORIGINAL) ──
        hc = HLA_CHAIN(chain_in, wrapper)
        chain = hc.matrices
        g3 = GATE_G3(chain.filter { a, _al, _re, _d -> a == 'original' }, wrapper)
        arm_mats = HLA_MATRICES(chain.filter { a, _al, _re, _d -> a != 'original' }
                                   .combine(channel.fromList(params.Cohorts))
                                   .combine(channel.fromList(params.HlaClasses)), src('hla_matrices.py')).matrix
        pub_mats = channel.fromList(params.Cohorts).combine(channel.fromList(params.HlaClasses)).map { c, k ->
            tuple('original', c, k, ext("${params.hla_v2}/results/06.assoc_inputs/${c}/${k}_dosage.tsv"))
        }

        // ── re-tests of every arm ────────────────────────────────────────────────
        gw = RETEST_GWAS(geno.combine(channel.of('best_guess', 'mean')), bedio, src('saige_collect.py'), lib).res
        gn_in = geno.combine(orig.map { c, _o, _gw, gn -> tuple(c, gn) }, by: 0)
            .map { c, a, gwf, gnf, ogn -> tuple(c, a, gwf, gnf, ogn) }
            .combine(channel.fromList(params.Strata))
        gn = RETEST_GENE(gn_in, genes, src('gene_retest.py'), bedio, lib)
        hla = RETEST_HLA(pub_mats.mix(arm_mats))
        omni_in = hla.plink.filter { _a, _c, k, _b, _m, _f -> k == 'residue' }
            .map { a, c, _k, b, m, f -> tuple(a, params.Cohorts.indexOf(c), c, b, m, f) }
            .groupTuple(size: params.Cohorts.size())
            .map { a, ix, cs, bs, ms, fs ->
                def o = (0..<ix.size()).toList().sort { i -> ix[i] }
                tuple(a, o.collect { i -> cs[i] }, o.collect { i -> bs[i] }, o.collect { i -> ms[i] }, o.collect { i -> fs[i] })
            }
        omni = RETEST_OMNIBUS(omni_in).res

        // ── gates G1-G4 on the ORIGINAL arm ──────────────────────────────────────
        g1 = GATE_G1(gw.filter { _c, a, imp, _t, _r -> a == 'original' && imp == 'best_guess' }, src('saige_collect.py')).json
        g4 = GATE_G4(hla.assoc.filter { a, _c, _k, _f -> a == 'original' }.map { _a, _c, _k, f -> f }.collect(),
                     omni.filter { a, _f -> a == 'original' }.map { _a, f -> f }, src('hla_gate.py'))
        GATE_SUMMARY(g1.mix(gn.gate, g3, g4).collect())

        // ── concordance of the arms (genotypes, QC-rule failures, HLA calls) ────
        if (doDs) {
            full = orig.filter { c, _a, _gw, _gn -> c == 'full_mainland' }.map { _c, _a, gwf, gnf -> tuple(gwf, gnf) }
            conc = CONCORDANCE(calls, full, targets.fixed, plan.plan, src('concordance.py'), bedio)
            QC_REASONS(calls, targets.fixed, plan.plan, table, src('concordance.py'), bedio)
            arm_specs = hc.evidence.filter { a, _p -> a != 'original' }.map { a, p -> "${a}=${p}" }.collect()
            mat_specs = arm_mats.map { a, c, k, p -> "${a}:${c}:${k}=${p}" }.collect()
            hcon = HLA_CONCORDANCE(arm_specs, mat_specs, plan.deep, table, src('hla_concordance.py'))
            GATES_DS(conc.summary, sub_depth.collect(), hcon.calls, src('gates_ds.py'))
        }

        // ── platform sensitivity: strata and leave-one-platform-out, shared controls ──
        if (doPf) {
            def subset_names = ['pooled'] + params.CasePlatformSlugs.collect { x -> "only_${x}" } +
                               params.CasePlatformSlugs.collect { x -> "without_${x}" }
            // control-source subsets: every case with one source's controls, or with all sources but one; they
            // run through the same null fits, step 2, gene tests and omnibus as the platform subsets
            def ctrl_names = params.ControlSources.collect { x -> "ctrlonly_${x}" } +
                             params.ControlSources.collect { x -> "ctrlwithout_${x}" }
            subs = MAKE_SUBSETS(channel.fromList(params.Cohorts).combine(table), src('make_subsets.py'), lib).dir
            csubs = MAKE_CTRL_SUBSETS(subs, ext(params.info_xlsx), src('control_subsets.py'), lib).dir
            null_in = subs.flatMap { c, d -> subset_names.collect { x -> tuple(c, x, file("${d}/${x}.random.tsv")) } }
                .mix(csubs.flatMap { c, d -> ctrl_names.collect { x -> tuple(c, x, file("${d}/${x}.random.tsv")) } })
            nm = FIT_NULL(null_in).null_model
            pg = PLATFORM_GENO(channel.fromList(params.Cohorts).combine(targets.fixed).combine(targets.null_snps), bedio, lib).geno
            hla_orig = hla.plink.filter { a, _c, _k, _b, _m, _f -> a == 'original' }
            h_al = hla_orig.filter { _a, _c, k, _b, _m, _f -> k == 'allele' }.map { _a, c, _k, b, m, f -> tuple(c, [b, m, f]) }
            h_re = hla_orig.filter { _a, _c, k, _b, _m, _f -> k == 'residue' }.map { _a, c, _k, b, m, f -> tuple(c, [b, m, f]) }
            s2 = STEP2_SUBSET(nm.combine(pg, by: 0).combine(h_al, by: 0).combine(h_re, by: 0),
                              bedio, src('saige_collect.py'), lib).res

            orig_gene = orig.map { c, _a, _gw, gnf -> tuple(c, gnf) }
            gene_in = subs.flatMap { c, d -> subset_names.collect { x -> tuple(c, x, file("${d}/${x}.rv_include.txt")) } }
                .mix(csubs.flatMap { c, d -> ctrl_names.collect { x -> tuple(c, x, file("${d}/${x}.rv_include.txt")) } })
                .combine(orig_gene, by: 0).combine(channel.fromList(params.Strata))
            PLATFORM_GENE(gene_in, genes, src('gene_retest.py'), bedio, lib)
            cf = CARRIER_FIRTH(orig_gene.join(subs).combine(genes), src('gene_retest.py'), src('carrier_firth.py'), bedio, lib).res

            om_base = h_re.map { c, files -> tuple(c, files[0], files[1], files[2]) }.join(subs).toList()
                .map { L -> L.sort { x -> params.Cohorts.indexOf(x[0]) } }
            om_base_c = h_re.map { c, files -> tuple(c, files[0], files[1], files[2]) }.join(csubs).toList()
                .map { L -> L.sort { x -> params.Cohorts.indexOf(x[0]) } }
            om_in = om_base.flatMap { L ->
                subset_names.collect { x ->
                    tuple(x, L.collect { y -> y[0] }, L.collect { y -> y[1] }, L.collect { y -> y[2] }, L.collect { y -> y[3] },
                          L.collect { y -> y[4] })
                }
            }.mix(om_base_c.flatMap { L ->
                ctrl_names.collect { x ->
                    tuple(x, L.collect { y -> y[0] }, L.collect { y -> y[1] }, L.collect { y -> y[2] }, L.collect { y -> y[3] },
                          L.collect { y -> y[4] })
                }
            })
            PLATFORM_OMNIBUS(om_in)

            stats_in = s2.filter { _c, x, _f -> x in subset_names }.groupTuple(size: subset_names.size()).join(subs).join(pg).join(cf)
            PLATFORM_STATS(stats_in, targets.null_snps, targets.signals, table, src('platform_stats.py'), bedio, lib)

            g5_in = nm.filter { _c, x, _r, _v -> x == 'pooled' }.map { c, _x, _r, v -> tuple(c, v) }
                .join(s2.filter { _c, x, _f -> x == 'pooled' }.map { c, _x, f -> tuple(c, f) })
                .join(gw.filter { _c, a, imp, _t, _r -> a == 'original' && imp == 'best_guess' }.map { c, _a, _i, t, _r -> tuple(c, t) })
            GATE_G5(g5_in)
        }
    }

    // ── report: summaries and figures from the published results ─────────────
    if ('report' in runs()) {
        // report_lib.py (shared helpers) and checks.py (the check registry) with the frozen modules they import
        def rlib = [lib, bedio, src('control_subsets.py'), src('report_lib.py'), src('checks.py')]
        t = DS_SUMMARY(src('ds_summary.py'), rlib).tables
        pr = PLATFORM_REPORT(src('platform_report.py'), rlib).tables
        // the per-cohort array comparison and HIBAG (heavy, cached), then the evidence per signal (seconds, rerun)
        acx = ARRAY_COHORTS([lib, src('array_cohorts.py'), src('array_hibag.py'), src('array_hibag.R')]).cohorts
        ae = ARRAY_EVIDENCE(src('array_evidence.py'), rlib, acx).tables
        cr = CONTROL_REPORT(src('control_report.py'), rlib + [src('platform_stats.py')]).tables
        HLA_REF_SHIFT(src('hla_reference_shift.py'), rlib)     // its tables are kept; no check reads them
        fu = FOLLOWUP(src('signal_followup.py'), rlib).tables
        ev = EVIDENCE(src('evidence_matrix.py'), rlib, t.mix(ae).mix(cr).mix(fu).collect()).tables
        FIGURES(t.mix(pr).mix(ae).mix(cr).mix(fu).mix(ev).collect(), src('figures.py'), src('catalogue.py'), rlib)
    }
}
