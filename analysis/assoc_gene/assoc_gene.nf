nextflow.enable.dsl = 2

// ═══════════════════════════════════════════════════════════════════════════
// assoc_gene — gene-based rare-variant association on one canonical variant
// map: rvtest CMC and SKAT-O under a fixed-effects logistic model, three nested
// cohorts, two nested impact strata, and a tiered robustness deliverable.
//
// Three design points, in the order they matter:
//
//   1. A GENE'S VARIANT SET IS A FILE. Gene membership comes from snpEff's own
//      gene attribution of each variant, materialised as map.<stratum>.tsv.gz
//      by MERGE_MAP, and rvtest is handed exactly that set through --setFile
//      (one 1-bp range per variant). So "which variants did this gene test"
//      is answerable from disk, per gene, which is what makes the per-gene
//      carrier counts and the publication summary table computable at all.
//   2. THE DENOMINATOR IS FIXED BEFORE ANY P-VALUE EXISTS. MERGE_MAP applies
//      the >= MinNumVar rule and counts the surviving genes; that count is the
//      Bonferroni denominator for every method in that cohort x stratum. BH
//      runs beside it over the identical family (genes rvtest did not return
//      enter at p = 1). Neither rule is corrected across methods or cohorts:
//      they test the same genes on the same samples, not independent
//      hypotheses. BH q < 0.05 is "called"; Bonferroni is the stricter tier.
//   3. ROBUSTNESS IS TIERED, NOT BINARY, AND IT IS NOT REPLICATION. The three
//      cohorts are strictly nested (full_mainland has 20 more cases than
//      narrow_mainland; the difference is controls), and the two strata are
//      nested. A gene called in all three cohorts has survived adding
//      controls. ROBUSTNESS counts that on the cohort axis, reports the tier,
//      and every figure sidecar says what the count means.
//
// Per cohort:
//
//   PREP_PHENO_COV ─────────────────────────────────────────────┐
//   SPLIT_BIM ─┐                                                │
//   MAKE_MINAC_KEEP ─┴→ MAP_CHUNK (x22) → MERGE_MAP  ← B1: the denominator
//        │                                    │
//        └───────────────→ EMIT_TEST_INPUTS ──┤  setFile + set/variant stats
//                                             ├→ RVTEST (x22 x 2 strata x 2 methods)
//                                             ↓
//                                       MERGE_RVTEST  ← B2
//                                             ↓
//                                           SCAN      ← B3: map identity + BH
//                                             ↓
//                              SIGNIFICANCE → PLOT_SCAN
//
// ...then across cohorts: COLLECT_ALL ← B4 → ROBUSTNESS → PLOT_GRID /
// PLOT_ROBUST_GENES / PLOT_GENE_DETAIL.
//
// Full methodology: docs/METHODS.md. Output dictionary: docs/OUTPUTS.md.
// Every check that gates the run: docs/VERIFICATION.md.
// ═══════════════════════════════════════════════════════════════════════════

// ── Cohorts ────────────────────────────────────────────────────────────────
params.Cohorts      = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
params.ValidCohorts = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']

// ── Paths ──────────────────────────────────────────────────────────────────
params.project_dir       = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6'
params.script_dir        = "${params.project_dir}/analysis/assoc_gene/scripts"
params.shared_script_dir = "${params.project_dir}/analysis/_shared/scripts"
params.out_dir           = "${params.project_dir}/analysis/assoc_gene/results"
params.model_inputs      = "${params.project_dir}/wgs.auto.par/results/12_model_inputs"
params.info_dir          = "${params.project_dir}/info"

// The genotypes are tuning.rv's fully-filtered call set, which IS
// 12_model_inputs/<c>/genotype/fixed_model/*.maf_lt_threshold after tuning.rv's
// sample --remove (3 QC outliers) and variant --exclude (the depth-diff knee).
// Verified a strict subset in both dimensions, so the QC provenance is
// tuning.rv's and is not repeated here.
params.TuningDir   = "${params.project_dir}/tuning.rv/results"
params.CallsetTpl  = "${params.TuningDir}/@@COHORT@@/02.callset_filter/filtered"
params.MinACTpl    = "${params.TuningDir}/@@COHORT@@/05.qc_collect/minac_recommendation.tsv"
params.MinAC       = null    // null = read tuning.rv's recommendation per cohort

params.PhenoTpl    = "${params.model_inputs}/@@COHORT@@/phenotype/pheno.tsv"
params.CovarTpl    = "${params.model_inputs}/@@COHORT@@/covariates/bbj_mainland_pc.sex.tsv"

// Pre-computed JHRP v6 snpEff annotations, one position-sorted tabix-indexed
// TSV per chromosome, keyed on (#CHROM, POS, REF, ALT). Our .bim IDs are
// literally chr<N>:<POS>:<REF>:<ALT>, so the join needs no conversion — which
// is why this component has no VCF annotation stage.
params.SnpEffIndex   = '/LARGE1/gr10478/platform/JHRPv6/workspace/pipeline/output/snpEff.v6.index'
params.SnpEffPattern = 'all.VQSR3.chr@@CHROM@@.vcf_out.tsv.2.gz'

// ── Tools ──────────────────────────────────────────────────────────────────
params.conda_env = 'cteph_geno_pro'
params.tool_dir  = '/home/b/b37974'
params.plink2    = '/home/b/b37974/plink2_alpha6/plink2'
params.bcftools  = '/home/b/b37974/bcftools/bcftools'
params.tabix     = '/home/b/b37974/htslib-1.9/tabix'
params.rvtestBin = '/home/b/b37974/rvtests/executable'

// ── The two variant sets ───────────────────────────────────────────────────
// MODIFIER is excluded from both: it is what snpEff gives intergenic and
// deep-intronic variants, which a gene-based test has no principled way to
// assign. Order is WIDEST FIRST: MAP_CHUNK annotates against the first entry's
// impact list and MERGE_MAP derives every narrower stratum from it, so the
// nesting is a property of the code path, not of the configuration.
params.Strata = [
    [ tag: 'low_moderate_high', impacts: 'HIGH MODERATE LOW' ],
    [ tag: 'moderate_high',     impacts: 'HIGH MODERATE'     ],
]

// ── The model ──────────────────────────────────────────────────────────────
// Logistic regression, sex + the first NPcs bbj_mainland principal components.
// rvtest reads --covar-name against the LOWER-CASED header, which is how
// prep_pheno_cov.py writes covar_rvt.tsv. Change NPcs (or CovarName outright)
// to change the model; RVTEST's pre-flight gate refuses a covariate the file
// does not carry, because rvtest itself would exit 0 on an uncovaried fit.
params.NPcs      = 10
params.PcLabel   = 'bbj_mainland'
params.PhenoName = 'PHENO1'
params.CovarName = 'sex,' + (1..params.NPcs).collect { "pc${it}_avg" }.join(',')

// Missing genotypes are imputed to the mean (2 x AF). It is rvtest's own
// default, passed explicitly so the log gate in RVTEST can assert it took
// effect: rvtest accepts an unrecognised --impute value SILENTLY and exits 0
// with unconsolidated genotypes (Main.cpp:1002-1013 has no else branch).
// For CMC the choice is immaterial -- the collapse casts genotypes to int, so
// an imputed 2*AF < 1 becomes 0 -- but the SKAT-O kernel uses the dosage.
params.ImputeMethod = 'mean'

// rvtest's two methods. `out` is the suffix rvtest actually puts on its result
// file: it writes <--out>.<TestName>.assoc, NOT <--out>.assoc, and the name is
// the test's own capitalisation rather than the flag's. Measured:
//     --burden cmc    -> chr4.moderate_high.cmc.CMC.assoc
//     --kernel skato  -> chr4.moderate_high.skato.SkatO.assoc
// Getting this wrong does not produce a wrong result, it produces no result:
// every RVTEST task dies on `mv: cannot stat`.
//
// `wald` is the second model rvtest runs in the SAME invocation, and it exists
// only to supply an effect size. --burden takes a comma-separated model list
// (ModelManager.cpp:37), so cmc and cmcWald share one tabix pass and one
// genotype extraction. cmcWald is a Wald fit of exactly the model the cmc
// SCORE test produces the published p from -- same cmcCollapse, same
// flipped-to-minor genotypes, same covariates -- so OR = exp(beta) is the
// model's own estimate, but z = beta/se does NOT reproduce the published p, and
// only the score p is published. See docs/METHODS.md.
//
// A KERNEL TEST ESTIMATES NO EFFECT SIZE, so skato carries wald: null and its
// RVTEST task writes an EMPTY chr<N>.wald.assoc rather than no file: an empty
// file keeps the DAG rectangular and lets merge_rvtest.py ASSERT emptiness per
// method instead of inferring it.
params.RvtestMethods = [
    [ tag: 'cmc',   opt: '--burden cmc,cmcWald', out: 'CMC',   wald: 'CMCWald' ],
    [ tag: 'skato', opt: '--kernel skato',       out: 'SkatO', wald: null      ],
]

// ── Map, significance, robustness ──────────────────────────────────────────
params.MinNumVar      = 2          // a gene needs at least this many mapped variants
params.Alpha          = 0.05
params.CollisionToken = '@'        // <symbol>@chr<N> when a symbol spans chromosomes
params.MaxGeneSpan    = 3000000
params.GeneAssignment = 'per_gene_most_severe'
params.StrictIdentity = true       // SCAN aborts on an unexplained set-size mismatch
params.MaxLabelGenes  = 24         // hard cap on the names in one Manhattan

// The gene model panel of the per-gene figures is drawn from Ensembl 86 through
// the shared region_tracks.R, which runs in the r_work environment.
params.rscript           = '/home/b/b37974/anaconda3/envs/r_work/bin/Rscript'

// The cohort whose numbers fill the publication summary table and the per-gene
// detail figures. All three cohorts' tables are written; this one is the one
// named in the figure captions. The largest sample is the natural choice.
params.ReportCohort          = 'full_mainland'
params.MaxRobustGenesPlotted = 40

params.Chroms = (1..22).collect { it.toString() }

// ═══════════════════════════════════════════════════════════════════════════
// Helpers
// ═══════════════════════════════════════════════════════════════════════════

def tpl(String t, String cohort)  { t.replace('@@COHORT@@', cohort) }
def snpeffFor(String chrom)       { "${params.SnpEffIndex}/" + params.SnpEffPattern.replace('@@CHROM@@', chrom) }
def script_file(String name)      { file("${params.script_dir}/${name}", checkIfExists: true) }

// Resolve the minAC at DAG-build time so the value that actually ran is knowable
// from the manifest rather than only from a log.
def minacFor(String cohort) {
    if (params.MinAC != null) return params.MinAC as Integer
    def f = file(tpl(params.MinACTpl, cohort))
    if (!f.exists()) {
        error "tuning.rv has not produced a minAC recommendation for '${cohort}'.\n" +
              "  expected: ${f}\n" +
              "  run tuning.rv first, or pass --MinAC <n> to override."
    }
    def row = f.readLines().find { ln -> ln.startsWith('recommended_minac') }
    if (!row) error "no 'recommended_minac' key in ${f}"
    return row.split('\t')[1].trim() as Integer
}

// Attach the chromosome parsed out of a fanned-out file name. Nextflow gives a
// task's multi-file output back as one list; the chromosome has to come from the
// name because the process cannot emit 22 separately-keyed tuples.
def byChrom(ch, String pattern) {
    ch.flatMap { c, files ->
        (files instanceof List ? files : [files]).collect { f ->
            def m = (f.name =~ pattern)
            if (!m) error "cannot parse a chromosome out of ${f.name} with ${pattern}"
            tuple(c, m[0][1].toString(), f)
        }
    }
}

// Pick one stratum's file out of a fanned-out list, by exact name.
def pick(files, String prefix, String stratum, String suffix) {
    def want = "${prefix}.${stratum}${suffix}"
    def hit = (files instanceof List ? files : [files]).find { f -> f.name == want }
    if (!hit) error "MERGE_MAP did not emit ${want}"
    return hit
}

// ═══════════════════════════════════════════════════════════════════════════
// STAGE 0 — inputs
// ═══════════════════════════════════════════════════════════════════════════

/* PREP_PHENO_COV — the phenotype and covariates on exactly the genotyped
 * samples, in the 1/2 coding rvtest requires, from one normalised vector.
 *
 * rvtest reads a 0 phenotype as MISSING (DataLoader.cpp:1248), so a 0/1 column
 * silently yields zero cases and a full result table with exit 0. The coding is
 * therefore normalised here, the counts are recorded in pheno_coding.json, and
 * RVTEST's post-flight gate checks the engine's own "Loaded N cases" line
 * against them. The case/control keep lists are written from the same vector so
 * the per-set allele counts and the association can never disagree about who is
 * a case. */
process PREP_PHENO_COV {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "${cohort}"

    publishDir { "${params.out_dir}/${cohort}/00.prep" }, mode: 'copy'

    input:
    tuple val(cohort), path(pheno), path(covar), path(fam)
    path script

    output:
    tuple val(cohort), path('pheno_rvt.tsv'), path('covar_rvt.tsv'),
          path('pheno_coding.json'), emit: pc
    tuple val(cohort), path('case.keep'), path('control.keep'), emit: keep

    script:
    """
    source activate ${params.conda_env}
    python3 ${script} \
        --pheno ${pheno} --covar ${covar} --fam ${fam} \
        --cohort ${cohort} --pheno-col ${params.PhenoName} \
        --out-rvt-pheno pheno_rvt.tsv \
        --out-rvt-covar covar_rvt.tsv \
        --out-case-keep case.keep --out-control-keep control.keep \
        --out-json pheno_coding.json
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// STAGE 1 — the map
// ═══════════════════════════════════════════════════════════════════════════

process SPLIT_BIM {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "${cohort}"

    input:
    tuple val(cohort), path(bim)

    output:
    tuple val(cohort), path('bim.chr*.tsv'), emit: chunks

    script:
    """
    set -euo pipefail
    awk -v OUT=. '\$1 ~ /^([1-9]|1[0-9]|2[0-2])\$/ { print > (OUT"/bim.chr"\$1".tsv") }' ${bim}
    # A chromosome with no variants would silently drop 1/22 of the genome.
    n=\$(ls bim.chr*.tsv | wc -l)
    test "\$n" -eq ${params.Chroms.size()} || { echo "ERROR: \$n autosome slices, expected ${params.Chroms.size()}" >&2; exit 1; }
    """
}

/* MAKE_MINAC_KEEP — the minAC decision as a variant-ID LIST: one artefact that
 * the map and the engine share, rather than two independent re-derivations. */
process MAKE_MINAC_KEEP {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag "${cohort}:mac${minac}"

    publishDir { "${params.out_dir}/${cohort}/01.map" }, mode: 'copy', pattern: '*.log'

    input:
    tuple val(cohort), val(minac), path(bed), path(bim), path(fam)

    output:
    tuple val(cohort), path('minac.snplist'), emit: keep
    path 'minac.log'

    script:
    """
    set -euo pipefail
    ${params.plink2} --bfile ${bed.baseName} --mac ${minac} \
        --write-snplist --out minac --threads ${task.ext.threads ?: 1}
    test -s minac.snplist
    """
}

/* MAP_CHUNK — one chromosome's variant -> gene map: a streaming sorted-merge
 * join of the .bim against the pre-computed JHRP snpEff index. Measured 14.4 s
 * for chr22. */
process MAP_CHUNK {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag "${cohort}:chr${chrom}"

    input:
    tuple val(cohort), val(chrom), path(bimchunk), path(keep)
    path script

    output:
    tuple val(cohort), path("map.chr${chrom}.tsv.gz"), emit: map
    tuple val(cohort), path("anno_counts.chr${chrom}.tsv"), emit: anno
    tuple val(cohort), path("map_qc.chr${chrom}.tsv"), emit: qc

    script:
    // --impacts must stay SHELL-UNQUOTED so argparse nargs='+' word-splits it.
    // Quoting it makes argparse receive one string, which matches nothing and
    // yields an EMPTY MAP — an empty scan, not an error. The widest stratum's
    // impact list is annotated; MERGE_MAP derives the narrower strata from it.
    def impacts = params.Strata[0].impacts
    """
    source activate ${params.conda_env}
    python3 ${script} \
        --bim ${bimchunk} --keep-variants ${keep} \
        --snpeff-index ${snpeffFor(chrom)} --chrom ${chrom} \
        --impacts ${impacts} \
        --gene-assignment ${params.GeneAssignment} \
        --out-map map.chr${chrom}.tsv.gz \
        --out-anno anno_counts.chr${chrom}.tsv \
        --out-qc map_qc.chr${chrom}.tsv
    """
}

/* MERGE_MAP — BARRIER B1. The denominator is fixed HERE.
 *
 * merge_gene_map.py applies the >= MinNumVar rule FIRST and then counts the
 * surviving genes, so n_genes_mapped -- and alpha / n_genes_mapped -- depends
 * only on the map and is known before a single association runs. The genes
 * dropped below the threshold are recorded (n_genes_dropped_lt_minvar) but
 * never counted. Both decision rules use this family: Bonferroni as its
 * threshold, BH as the size of the family it corrects over. */
process MERGE_MAP {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag "${cohort}"

    publishDir { "${params.out_dir}/${cohort}/01.map" }, mode: 'copy'

    input:
    tuple val(cohort), path(maps, stageAs: 'map_*.tsv.gz'),
          path(annos, stageAs: 'anno_*.tsv'), path(qcs, stageAs: 'qc_*.tsv')
    path script

    output:
    tuple val(cohort), path('map.*.tsv.gz'), path('gene_index.*.tsv'),
          path('variants.*.txt'), path('denominator.*.tsv'),
          path('map_collisions.*.tsv'), emit: strata
    tuple val(cohort), path('anno_counts.tsv'), path('map_qc.tsv'), emit: anno

    script:
    def strata_args = params.Strata.collect { st -> "--stratum ${st.tag}=${st.impacts.replaceAll(' ', ',')}" }.join(' ')
    """
    source activate ${params.conda_env}
    python3 ${script} \
        --map ${maps} --anno ${annos} --qc ${qcs} \
        ${strata_args} \
        --min-num-var ${params.MinNumVar} --alpha ${params.Alpha} \
        --collision-token '${params.CollisionToken}' \
        --max-gene-span ${params.MaxGeneSpan} \
        --cohort ${cohort} --out-dir .
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// STAGE 2 — test inputs
// ═══════════════════════════════════════════════════════════════════════════

/* EMIT_TEST_INPUTS — the tested genotypes, the rvtest set definitions and the
 * per-gene allele and carrier counts, in ONE task from ONE map.
 *
 * The plink2 calls and the file emission are deliberately not separate
 * processes. The chromosome vocabulary written into the setFile is READ BACK
 * from the VCF that was just written, and the variant IDs are read back from
 * the .bim as plink2 wrote it — so a stale bfile can never be paired with a
 * fresh set definition, because there is no channel between them.
 *
 * THE PER-GENE COUNTS ARE COMPUTED HERE FOR THE SAME REASON. rvtest reports no
 * allele count and no carrier count for a set test. Both are grouped with the
 * SAME in-memory variant_id -> set_name dict that just wrote the set file:
 * cumulative MAC from two plink2 --freq counts runs (one per group), carriers
 * from a direct decode of tested.bed -- and the bed-derived MAC must equal the
 * .acount-derived MAC for EVERY variant, which is the cross-check that proves
 * the decode read the right allele in the right sample order. */
process EMIT_TEST_INPUTS {
    executor 'slurm'
    queue 'gr10478b'
    time '8h'
    tag "${cohort}:${stratum}"

    publishDir { "${params.out_dir}/${cohort}/02.test_inputs/${stratum}" }, mode: 'copy',
               pattern: '{rvtest.set.chr*.txt,test_inputs_qc.tsv,set_stats.tsv,variant_stats.tsv,multi_carriers.tsv}'

    input:
    tuple val(cohort), val(stratum), path(map), path(gene_index), path(variants),
          path(bed), path(bim), path(fam), path(case_keep), path(control_keep),
          path(coding)
    path script

    output:
    tuple val(cohort), val(stratum), path('tested.vcf.gz'), path('tested.vcf.gz.tbi'), emit: vcf
    tuple val(cohort), val(stratum), path('rvtest.set.chr*.txt'), emit: setfiles
    tuple val(cohort), val(stratum), path('set_stats.tsv'), emit: setstats
    tuple val(cohort), val(stratum), path('variant_stats.tsv'), emit: varstats
    tuple val(cohort), val(stratum), path('multi_carriers.tsv'), emit: multi
    tuple val(cohort), val(stratum), path('test_inputs_qc.tsv'), emit: qc

    script:
    """
    export PATH=${params.tool_dir}:\$PATH
    source activate ${params.conda_env}
    set -euo pipefail

    ${params.plink2} --bfile ${bed.baseName} --extract ${variants} \
        --make-bed --out tested --threads ${task.ext.threads ?: 1}
    ${params.plink2} --bfile tested --export vcf id-paste=iid bgz \
        --out tested --threads ${task.ext.threads ?: 1}
    ${params.tabix} -f -p vcf tested.vcf.gz

    # Every sample must survive; only variants are being selected.
    test \$(wc -l < tested.fam) -eq \$(wc -l < ${fam})

    # --- per-set allele counts, on exactly the TESTED samples ----------------
    # --freq counts is founders-only by default. Every sample here is a founder,
    # so --nonfounders is a no-op today -- and it makes the COUNTED sample set
    # equal the TESTED one by construction rather than by a property of the .fam
    # that upstream could change without this noticing.
    test \$(awk '\$3!="0" || \$4!="0"' tested.fam | wc -l) -eq 0
    # The two groups PARTITION the sample set (PREP_PHENO_COV aborts on any
    # missing phenotype), so the pooled counts are the sums and a third plink2
    # call would be a second chance to disagree with them.
    test \$(( \$(wc -l < ${case_keep}) + \$(wc -l < ${control_keep}) )) -eq \$(wc -l < tested.fam)

    ${params.plink2} --bfile tested --keep ${case_keep}    --nonfounders --freq counts \
        --out mac_case    --threads ${task.ext.threads ?: 1}
    ${params.plink2} --bfile tested --keep ${control_keep} --nonfounders --freq counts \
        --out mac_control --threads ${task.ext.threads ?: 1}

    python3 ${script} \
        --map ${map} --gene-index ${gene_index} \
        --bed tested.bed --bim tested.bim --fam tested.fam --vcf tested.vcf.gz \
        --bcftools ${params.bcftools} \
        --acount-case mac_case.acount --acount-control mac_control.acount \
        --case-keep ${case_keep} --control-keep ${control_keep} --coding ${coding} \
        --cohort ${cohort} --stratum ${stratum} \
        --out-dir . --out-qc test_inputs_qc.tsv \
        --out-set-stats set_stats.tsv --out-variant-stats variant_stats.tsv \
        --out-multi-carriers multi_carriers.tsv
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// STAGE 3 — the associations
// ═══════════════════════════════════════════════════════════════════════════

/* RVTEST — one task per cohort x stratum x method x chromosome.
 *
 * --setFile, never --geneFile. rvtest errors if both are given (Main.cpp:919),
 * which is the guard that no refFlat-span grouping is in play. NOTE that
 * --setFile also renames the first output column from 'Gene' to 'Range'; that is
 * absorbed in MERGE_RVTEST, which asserts the header rather than trusting it.
 *
 * TWO GATES, both non-optional:
 *   pre-flight  — rvtest EXITS 0 when a --covar-name column is absent and
 *                 silently fits an UNCOVARIED model.
 *   post-flight — rvtest EXITS 0 when the phenotype is 0/1-coded and silently
 *                 analyses zero cases. Measured: "Loaded 0 cases, 418 controls,
 *                 1747 missing", a full result table, exit 0. A downstream check
 *                 cannot tell that from a genuinely null scan, so it is gated
 *                 HERE against the counts PREP_PHENO_COV recorded.  */
process RVTEST {
    executor 'slurm'
    queue 'gr10478b'
    time '8h'
    tag "${cohort}:${stratum}:${method}:chr${chrom}"

    publishDir { "${params.out_dir}/${cohort}/03.assoc/${stratum}/${method}" }, mode: 'copy',
               pattern: '*.rvtest.log'

    input:
    tuple val(cohort), val(stratum), val(method), val(opt), val(suffix), val(wald),
          val(chrom), path(vcf), path(tbi), path(setfile), path(pheno), path(covar),
          path(coding)

    output:
    tuple val(cohort), val(stratum), val(method), val(chrom),
          path("chr${chrom}.assoc"), path("chr${chrom}.wald.assoc"), emit: assoc
    path "chr${chrom}.rvtest.log"

    script:
    """
    export PATH=${params.rvtestBin}:\$PATH
    set -euo pipefail

    # --- pre-flight: every named covariate must exist -----------------------
    head -1 ${covar} | tr '\\t ' '\\n' > .covar_header
    for c in \$(echo '${params.CovarName}' | tr ',' ' '); do
        grep -qx "\$c" .covar_header || {
            echo "ERROR: covariate column '\$c' is absent from ${covar}." >&2
            echo "rvtest would exit 0 and fit an uncovaried model." >&2
            exit 1
        }
    done

    rvtest \
        --inVcf ${vcf} \
        --pheno ${pheno} --pheno-name ${params.PhenoName.toLowerCase()} \
        --covar ${covar} --covar-name ${params.CovarName} \
        --setFile ${setfile} \
        --out chr${chrom}.${stratum}.${method} \
        --impute ${params.ImputeMethod} \
        --noweb --numThread ${task.ext.threads ?: 1} \
        ${opt} 2>&1 | tee chr${chrom}.rvtest.log

    # --- post-flight: the engine's own count of what it analysed ------------
    # grep -a is NOT optional: rvtest's log carries bytes grep classifies as
    # binary, so a plain grep returns "Binary file ... matches" INSTEAD of the
    # line and every numeric test below fails on a non-integer. Measured.
    exp_case=\$(python3 -c "import json;print(json.load(open('${coding}'))['n_case'])")
    exp_ctrl=\$(python3 -c "import json;print(json.load(open('${coding}'))['n_control'])")
    exp_n=\$(python3 -c "import json;print(json.load(open('${coding}'))['n_samples'])")

    if grep -aqi 'There are no case' chr${chrom}.rvtest.log; then
        echo "ERROR: rvtest analysed ZERO cases. This is the 0/1 phenotype defect." >&2
        exit 1
    fi
    line=\$(grep -a -m1 -E 'Loaded [0-9]+ cases' chr${chrom}.rvtest.log || true)
    test -n "\$line" || { echo "ERROR: rvtest never reported a case/control count." >&2; exit 1; }
    got_case=\$(echo "\$line" | sed -E 's/.*Loaded ([0-9]+) cases.*/\\1/')
    got_ctrl=\$(echo "\$line" | sed -E 's/.*, ([0-9]+) controls.*/\\1/')
    got_miss=\$(echo "\$line" | sed -E 's/.*and ([0-9]+) missing.*/\\1/')
    test "\$got_case" -eq "\$exp_case" || { echo "ERROR: \$got_case cases, expected \$exp_case" >&2; exit 1; }
    test "\$got_ctrl" -eq "\$exp_ctrl" || { echo "ERROR: \$got_ctrl controls, expected \$exp_ctrl" >&2; exit 1; }
    test "\$got_miss" -eq 0            || { echo "ERROR: \$got_miss missing phenotypes" >&2; exit 1; }

    # --- post-flight: --impute mean actually reached the consolidator --------
    # Main.cpp:1002-1013 matches the value against four branches and HAS NO ELSE.
    # An unrecognised value (--impute Mean, say) leaves DataConsolidator::strategy
    # at UNINITIALIZED (DataConsolidator.cpp:145), whose consolidate() branch only
    # calls logger->error(...) -- rvtest then finishes with EXIT 0 and
    # unconsolidated genotypes. There is no way to detect that downstream, so it
    # is asserted from the log here.
    if grep -aq 'Impute missing genotype to mean (by default)' chr${chrom}.rvtest.log; then
        echo "ERROR: rvtest fell back to its default; --impute never reached the parser." >&2
        exit 1
    fi
    grep -aq 'Impute missing genotype to mean' chr${chrom}.rvtest.log || {
        echo "ERROR: rvtest reported no imputation strategy at all." >&2
        exit 1
    }
    if grep -aq 'Uninitialized consolidation methods' chr${chrom}.rvtest.log; then
        echo "ERROR: --impute value was not recognised. rvtest logs this PER GENE and" >&2
        echo "still exits 0, leaving genotypes unconsolidated." >&2
        exit 1
    fi

    # Kernel weights are the engine default and are NOT passed on the command
    # line; this asserts the default is what we think it is. rvtest prints it
    # unconditionally for skato (ModelManager.cpp:186).
    if [ "${method}" = "skato" ]; then
        grep -aq 'beta1 = 1.00, beta2 = 25.00' chr${chrom}.rvtest.log || {
            echo "ERROR: SKAT-O did not report Beta(1,25) weights." >&2
            grep -a 'weight' chr${chrom}.rvtest.log >&2 || true
            exit 1
        }
    fi

    # rvtest names its output <--out>.<TestName>.assoc; see params.RvtestMethods.
    produced=chr${chrom}.${stratum}.${method}.${suffix}.assoc
    test -f "\$produced" || {
        echo "ERROR: rvtest did not write \$produced." >&2
        echo "It wrote:" >&2; ls -1 *.assoc 2>/dev/null >&2 || echo "  (no .assoc at all)" >&2
        exit 1
    }
    mv "\$produced" chr${chrom}.assoc

    # The second model of the same invocation, when this method has one. It
    # supplies beta/se only; MERGE_RVTEST asserts the block height and takes the
    # first row of each Range block (Model.h:961 writes 1 + n_covar rows per gene
    # over X = [Intercept | collapsed | covariates], and i = 1 is the collapsed
    # genotype). An empty file for a kernel method is the deliberate signal that
    # this method estimates no effect size -- see params.RvtestMethods.
    if [ -n "${wald}" ]; then
        w=chr${chrom}.${stratum}.${method}.${wald}.assoc
        test -s "\$w" || {
            echo "ERROR: rvtest did not write \$w." >&2
            echo "It wrote:" >&2; ls -1 *.assoc 2>/dev/null >&2 || echo "  (none)" >&2
            exit 1
        }
        mv "\$w" chr${chrom}.wald.assoc
    else
        : > chr${chrom}.wald.assoc
    fi

    # N_INFORMATIVE on the first tested gene must be the whole sample set.
    n_inf=\$(awk -F'\\t' 'NR==2{print \$3}' chr${chrom}.assoc)
    test "\$n_inf" -eq "\$exp_n" || { echo "ERROR: N_INFORMATIVE \$n_inf != \$exp_n samples" >&2; exit 1; }
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// STAGE 4 — merge, scan, significance
// ═══════════════════════════════════════════════════════════════════════════

/* MERGE_RVTEST — BARRIER B2. The 22 chromosome parts of one method become one
 * table: Range -> set_name, rho kept from SKAT-O, beta/se from cmcWald with the
 * OR and its 95 % CI derived once here, and the per-gene counts joined from
 * set_stats.tsv. */
process MERGE_RVTEST {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "${cohort}:${stratum}:${method}"

    publishDir { "${params.out_dir}/${cohort}/03.assoc/${stratum}/${method}" }, mode: 'copy'

    input:
    tuple val(cohort), val(stratum), val(method), val(chroms),
          path(parts, stageAs: 'part_*.assoc'),
          path(walds, stageAs: 'wald_*.assoc'), path(gene_index), path(set_stats)
    path script

    output:
    tuple val(cohort), val(stratum), val(method), path('merged.tsv'), emit: merged

    script:
    """
    source activate ${params.conda_env}
    python3 ${script} \
        --assoc ${parts} --wald ${walds} --gene-index ${gene_index} \
        --set-stats ${set_stats} \
        --n-covar ${params.CovarName.split(',').size()} \
        --cohort ${cohort} --stratum ${stratum} --method ${method} \
        --out merged.tsv
    """
}

/* SCAN — BARRIER B3. Both methods' tables become one scan per cohort x stratum.
 *
 * Two things happen here and nowhere else:
 *   1. THE MAP IDENTITY. For every gene, rvtest's NumVar must equal the map's
 *      n_var_map -- or exceed it by exactly the number of multi-base-REF records
 *      whose span covers a later mapped variant (a 1-bp range matches such a
 *      record twice; scan.py predicts the excess exactly, 81,746/81,746 genes
 *      on the reference data). Under --strict anything else aborts.
 *   2. THE TWO DECISION RULES. Bonferroni at the threshold MERGE_MAP wrote, and
 *      BH over the identical n_genes_mapped family with missing genes at p = 1.
 *      Because they share alpha and n, every Bonferroni call is a BH call;
 *      scan.py asserts it. */
process SCAN {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "${cohort}:${stratum}"

    publishDir { "${params.out_dir}/${cohort}/04.scan/${stratum}" }, mode: 'copy'

    input:
    tuple val(cohort), val(stratum), path(merged, stageAs: 'merged_*.tsv'),
          path(gene_index), path(denominator), path(collisions), path(map)
    path script

    output:
    tuple val(cohort), val(stratum), path('gene_scan.tsv'), path('scan_qc.tsv'), emit: scan
    tuple val(cohort), val(stratum), path('map_identity.tsv'), emit: identity

    script:
    def strict = params.StrictIdentity ? '--strict' : '--no-strict'
    """
    source activate ${params.conda_env}
    python3 ${script} \
        --assoc ${merged} \
        --gene-index ${gene_index} --denominator ${denominator} \
        --collisions ${collisions} --map ${map} \
        --cohort ${cohort} --stratum ${stratum} ${strict} \
        --out-scan gene_scan.tsv --out-qc scan_qc.tsv \
        --out-identity map_identity.tsv
    """
}

/* SIGNIFICANCE — one table per cohort x stratum x method stating the decision
 * rules that were actually applied: Bonferroni (recomputed, then checked
 * against the threshold the scan was judged at), BH, and nominal for reference. */
process SIGNIFICANCE {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "${cohort}:${stratum}:${method}"

    publishDir { "${params.out_dir}/${cohort}/05.signals" }, mode: 'copy',
               saveAs: { _f -> "significance.${method}.${stratum}.tsv" }

    input:
    tuple val(cohort), val(stratum), val(method), path(scan_qc), path(denominator)
    path script

    output:
    tuple val(cohort), val(stratum), val(method), path('significance.tsv'), emit: sig

    script:
    """
    source activate ${params.conda_env}
    python3 ${script} \
        --scan-qc ${scan_qc} --denominator ${denominator} \
        --cohort ${cohort} --stratum ${stratum} --method ${method} \
        --alpha ${params.Alpha} \
        --out significance.tsv
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// STAGE 5 — figures, per cohort
// ═══════════════════════════════════════════════════════════════════════════

// One y axis per variant set, decided once. The three cohorts of a set are the
// same statistic on nested data, so a reader compares point heights across their
// three figures -- which only works if the three share a scale, and a figure
// drawn from its own cell cannot know the other two. Only the DATA maximum is
// shared: each figure still draws its own threshold lines, whose denominators
// differ.
process SCAN_SCALE {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'scan_scale'

    publishDir "${params.out_dir}/_comparison/tables", mode: 'copy'

    input:
    path scans, stageAs: 'scan_*.tsv'
    path script

    output:
    path 'scan_scale.tsv', emit: scale

    script:
    """
    source activate ${params.conda_env}
    python3 ${script} --gene-scan ${scans} --out-tsv scan_scale.tsv
    """
}


process PLOT_SCAN {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "${cohort}:${stratum}"

    publishDir "${params.out_dir}/figures/02.gene_scan", mode: 'copy', pattern: '*.png'

    input:
    tuple val(cohort), val(stratum), path(scan), path(scan_qc), path(denominator)
    path scale
    path script

    output:
    path "gene_scan.${cohort}.${stratum}.png"
    // The member's numbers; CATALOGUE turns the six of them into the family README.
    path "gene_scan.${cohort}.${stratum}.stats.json", emit: stats

    script:
    """
    source activate ${params.conda_env}
    python3 ${script} \
        --gene-scan ${scan} --scan-qc ${scan_qc} --denominator ${denominator} \
        --cohort ${cohort} --stratum ${stratum} \
        --max-label-genes ${params.MaxLabelGenes} \
        --scan-scale ${scale} \
        --out-png gene_scan.${cohort}.${stratum}.png
    """
}

/* PLOT_ANNOTATION — the variant-annotation distribution per cohort. Not part of
 * the association and never a decision input: it is the picture of what the
 * call set is made of, and of how much of it each stratum keeps. */
process PLOT_ANNOTATION {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "${cohort}"

    publishDir "${params.out_dir}/figures/01.variant_sets", mode: 'copy',
               pattern: '{*.png,tables/**}'

    input:
    // NO stageAs here, deliberately. plot_annotation.py derives the stratum from
    // the FILE NAME -- gene_index.<stratum>.tsv carries no stratum column -- so
    // renaming them would destroy the only thing that pairs a gene index with
    // its own denominator.
    tuple val(cohort), path(anno), path(map_qc), path(gene_index), path(denom)
    path script

    output:
    path "variant_sets.${cohort}.png"
    path "variant_sets.${cohort}.stats.json", emit: stats
    path "tables/${cohort}/*.tsv"

    script:
    """
    source activate ${params.conda_env}
    python3 ${script} \
        --anno-counts ${anno} --map-qc ${map_qc} \
        --gene-index ${gene_index} --denominator ${denom} \
        --cohort ${cohort} \
        --out-png variant_sets.${cohort}.png --out-dir tables/${cohort}
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// STAGE 6 — across cohorts
// ═══════════════════════════════════════════════════════════════════════════

/* COLLECT_ALL — BARRIER B4. Every table that already carries its own
 * cohort/stratum key is concatenated as-is; cohort_samples.tsv is built from
 * the three pheno_coding.json files so the summary table's denominators come
 * from the same contract the RVTEST gate checked the engine against. */
process COLLECT_ALL {
    executor 'local'
    tag 'all'

    publishDir "${params.out_dir}/_comparison/tables", mode: 'copy'

    input:
    path scans,    stageAs: 'scan_*.tsv'
    path qcs,      stageAs: 'qc_*.tsv'
    path sigs,     stageAs: 'sig_*.tsv'
    path idents,   stageAs: 'ident_*.tsv'
    path denoms,   stageAs: 'den_*.tsv'
    path setstats, stageAs: 'setstats_*.tsv'
    path varstats, stageAs: 'varstats_*.tsv'
    path multis,   stageAs: 'multi_*.tsv'
    path codings,  stageAs: 'coding_*.json'

    output:
    tuple path('gene_scan_all.tsv'), path('scan_qc_all.tsv'),
          path('significance_all.tsv'), path('map_identity_all.tsv'),
          path('denominators.tsv'), path('set_stats_all.tsv'),
          path('variant_stats_all.tsv'), path('multi_carriers_all.tsv'),
          path('cohort_samples.tsv'), emit: tables

    script:
    """
    set -euo pipefail
    cat_one () {  # header from the first file only
        out=\$1; shift
        first=1
        for f in "\$@"; do
            if [ \$first -eq 1 ]; then head -1 "\$f" > "\$out"; first=0; fi
            tail -n +2 "\$f" >> "\$out"
        done
    }
    cat_one gene_scan_all.tsv     ${scans}
    cat_one scan_qc_all.tsv       ${qcs}
    cat_one significance_all.tsv  ${sigs}
    cat_one map_identity_all.tsv  ${idents}
    cat_one denominators.tsv      ${denoms}
    cat_one set_stats_all.tsv     ${setstats}
    cat_one variant_stats_all.tsv ${varstats}
    cat_one multi_carriers_all.tsv ${multis}

    python3 - ${codings} <<'PY'
import json, sys
rows = []
for p in sys.argv[1:]:
    d = json.load(open(p))
    rows.append((d['cohort'], d['n_samples'], d['n_case'], d['n_control']))
rows.sort()
with open('cohort_samples.tsv', 'w') as fh:
    fh.write('cohort\\tn_samples\\tn_case\\tn_control\\n')
    for r in rows:
        fh.write('\\t'.join(str(x) for x in r) + '\\n')
PY
    """
}

/* ROBUSTNESS — the deliverable. Every gene x cohort x stratum x method in one
 * evidence table; a tier per (gene, stratum); and the publication summary table
 * per cohort. It introduces NO threshold and NO p-value: the scan already paid
 * alpha / n_genes_mapped, and this asks the different question "does the reading
 * survive a change of statistic and of sample selection". See
 * scripts/robust_genes.py for the tier rule and for why the count is taken on
 * the cohort axis while the strata are carried as sensitivity columns. */
process ROBUSTNESS {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'robustness'

    publishDir "${params.out_dir}/_comparison/robust_genes", mode: 'copy'

    input:
    // cohort_samples.tsv is re-emitted under its own name so PLOT_GENE_DETAIL
    // gets it from one tuple; the input must therefore be staged under ANOTHER
    // name, or `cp` refuses to copy a file onto itself.
    tuple path(scan_all), path(qc_all), path(sig_all), path(ident_all), path(denoms),
          path(setstats_all), path(varstats_all), path(multi_all),
          path(cohort_samples, stageAs: 'cohort_samples.in.tsv')
    path gene_index_all
    path script

    output:
    tuple path('evidence.tsv'), path('tiers.tsv'), path('variant_detail.tsv'),
          path('multi_carrier_detail.tsv'), path('cohort_samples.tsv'), emit: robust
    path 'gene_windows.tsv', emit: windows
    path 'summary_table.*.tsv'
    path 'summary_table.md'

    script:
    def cohorts = params.Cohorts instanceof List ? params.Cohorts.join(' ') : params.Cohorts
    def strata  = params.Strata.collect { st -> st.tag }.join(' ')
    """
    source activate ${params.conda_env}
    python3 ${script} \
        --gene-scan-all ${scan_all} --scan-qc-all ${qc_all} \
        --gene-index-all ${gene_index_all} --cohort-samples ${cohort_samples} \
        --denominators ${denoms} --set-stats-all ${setstats_all} \
        --variant-stats-all ${varstats_all} --multi-carriers-all ${multi_all} \
        --cohorts ${cohorts} --strata ${strata} \
        --report-cohort ${params.ReportCohort} --alpha ${params.Alpha} \
        --out-evidence evidence.tsv --out-tiers tiers.tsv \
        --out-variant-detail variant_detail.tsv \
        --out-multi-detail multi_carrier_detail.tsv \
        --out-gene-windows gene_windows.tsv \
        --out-summary-prefix summary_table
    cp ${cohort_samples} cohort_samples.tsv
    """
}

process PLOT_GRID {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag 'grid'

    publishDir "${params.out_dir}/figures/03.calibration", mode: 'copy'

    input:
    tuple path(scan_all), path(qc_all)
    path script

    output:
    path 'calibration.png'
    path 'README.md'

    script:
    """
    source activate ${params.conda_env}
    python3 ${script} --gene-scan-all ${scan_all} --scan-qc-all ${qc_all} \
        --out-png calibration.png --out-md README.md
    """
}

process PLOT_ROBUST_GENES {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag 'robust_genes'

    publishDir "${params.out_dir}/figures/04.robust_genes", mode: 'copy'

    input:
    tuple path(evidence), path(tiers), path(cohort_samples), path(denoms)
    path script

    output:
    path 'robust_genes.png'
    path 'README.md'

    script:
    """
    source activate ${params.conda_env}
    python3 ${script} --evidence ${evidence} --tiers ${tiers} \
        --cohort-samples ${cohort_samples} --denominators ${denoms} \
        --report-cohort ${params.ReportCohort} --max-genes ${params.MaxRobustGenesPlotted} \
        --alpha ${params.Alpha} \
        --out-png robust_genes.png --out-md README.md
    """
}

/* GENE_MODELS — the exon/intron structure of every Tier 1/2 gene, from the
 * shared region_tracks.R (Ensembl 86, GRCh38), one window per gene as
 * ROBUSTNESS computed it. One task, not one per gene: the count is
 * data-dependent and the DAG's cardinality is asserted by verify.sh.
 * A gene the annotation has no model for leaves an EMPTY file and the figure
 * says so, so one missing symbol cannot cost the whole family its figures. */
process GENE_MODELS {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'gene_models'

    input:
    path windows
    path region_script
    path gene_helper          // region_tracks.R sources it from its own directory

    output:
    path 'gene_models', emit: models

    script:
    """
    set -euo pipefail
    mkdir -p gene_models
    tail -n +2 ${windows} | while IFS=\$'\\t' read -r gene set_name chrom start end; do
        [ -n "\${gene:-}" ] || continue
        out=gene_models/"\${gene}".exons.tsv
        ${params.rscript} ${region_script} --chrom "\${chrom}" --start "\${start}" \
            --end "\${end}" --out-exons "\${out}" --out-recomb /dev/null || :
        [ -s "\${out}" ] || : > "\${out}"
    done
    n_all=\$(ls gene_models | wc -l)
    n_ok=\$(find gene_models -name '*.exons.tsv' -size +0 | wc -l)
    echo "[gene_models] \${n_ok}/\${n_all} gene(s) with an Ensembl 86 model"
    if [ "\${n_all}" -gt 0 ] && [ "\${n_ok}" -eq 0 ]; then
        echo "[gene_models] ERROR: no gene got a model - R environment or staging is broken" >&2
        exit 1
    fi
    """
}

/* CATALOGUE — one README per fan-out figure family. The members drop their
 * numbers as <name>.stats.json; the prose is written once, here, with every
 * member as a row of one table. */
process CATALOGUE {
    executor 'local'
    tag "${family}"

    publishDir { "${params.out_dir}/figures/${dir}" }, mode: 'copy'

    input:
    tuple val(family), val(dir), path(stats, stageAs: 'stats_*.json')
    path script

    output:
    path 'README.md'

    script:
    """
    source activate ${params.conda_env}
    python3 ${script} --family ${family} --stats ${stats} --figure-dir ${dir} \
        --alpha ${params.Alpha} --min-num-var ${params.MinNumVar} --out README.md
    """
}

/* PLOT_GENE_DETAIL — one figure per Tier 1/2 gene. The PNGs are optional
 * outputs (a null scan tiers nothing); the index is mandatory so the DAG stays
 * rectangular and a run with no figures still says so on disk. */
process PLOT_GENE_DETAIL {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag 'gene_detail'

    publishDir "${params.out_dir}/figures/05.genes", mode: 'copy'

    input:
    tuple path(tiers), path(evidence), path(variant_detail), path(multi_detail),
          path(cohort_samples)
    path gene_models
    path script

    output:
    path '*.png', optional: true
    path 'README.md'
    path 'index.tsv'

    script:
    """
    source activate ${params.conda_env}
    python3 ${script} --tiers ${tiers} --evidence ${evidence} \
        --variant-detail ${variant_detail} --multi-detail ${multi_detail} \
        --cohort-samples ${cohort_samples} \
        --report-cohort ${params.ReportCohort} --exon-dir ${gene_models} \
        --out-dir . --out-index index.tsv --out-md README.md
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
// EVERY `by:` below is a plain-String index list, never a composite tuple key.
// A tuple used as a `by:` element matches nothing and the channel simply never
// emits — no error, no warning, and every downstream process is skipped.
//
// That hazard is real, but an .ifEmpty { error ... } clause is NOT the guard
// for it: Nextflow attaches ifEmpty's closure to the operator's afterStop,
// which also runs when Session.abort() force-terminates the operators, so it
// fires on every abort and replaces the real error in the console. The
// empty-channel hazard is instead covered post-hoc: results/_run_info/trace.txt
// carries one row per task, so a barrier that emitted nothing is a MISSING
// PROCESS NAME, and verify.sh checks the per-process task counts against this
// DAG's cardinality.
//
// combine(channel.fromList(<list of lists>)) SPREADS each inner list into the
// tuple; every destructuring closure below is written for the spread form.
// ═══════════════════════════════════════════════════════════════════════════

workflow {

    // ── validation, before anything is submitted ───────────────────────────
    def cohortList = (params.Cohorts instanceof List)
        ? params.Cohorts
        : params.Cohorts.toString().split(',').collect { x -> x.trim() }.findAll { x -> x }
    def bad = cohortList.findAll { ch -> !(ch in params.ValidCohorts) }
    if (bad) error "unknown cohort(s): ${bad.join(', ')}. Valid: ${params.ValidCohorts.join(', ')}"
    if (!(params.ReportCohort in cohortList)) {
        error "ReportCohort '${params.ReportCohort}' is not among the cohorts being run: ${cohortList}"
    }

    // RVTEST's post-flight gate greps rvtest's log for the literal
    // "Impute missing genotype to mean". rvtest words the line differently for
    // every strategy, so the gate cannot be templated -- pin the parameter here
    // instead so the two cannot drift apart silently.
    if (params.ImputeMethod != 'mean') {
        error "ImputeMethod is '${params.ImputeMethod}'. Only 'mean' is supported: the RVTEST\n" +
              "  log gate is written against rvtest's exact wording for it, and rvtest accepts\n" +
              "  an unrecognised --impute value SILENTLY (Main.cpp:1002-1013 has no else) and\n" +
              "  exits 0 with unconsolidated genotypes, so this must not be a free-form string."
    }

    def minac = cohortList.collectEntries { ch -> [(ch): minacFor(ch)] }
    log.info "assoc_gene | cohorts: ${cohortList.join(', ')}"
    log.info "assoc_gene | minAC (tuning.rv): ${minac}"
    log.info "assoc_gene | strata: ${params.Strata.collect { st -> st.tag }.join(', ')}"
    log.info "assoc_gene | model: logistic, covariates ${params.CovarName}"
    log.info "assoc_gene | min variants per gene: ${params.MinNumVar} | alpha: ${params.Alpha}"

    ch_cohort = channel.fromList(cohortList)

    ch_callset = ch_cohort.map { c ->
        def p = tpl(params.CallsetTpl, c)
        tuple(c, file("${p}.bed", checkIfExists: true),
                 file("${p}.bim", checkIfExists: true),
                 file("${p}.fam", checkIfExists: true))
    }

    // ── STAGE 0 ────────────────────────────────────────────────────────────
    ch_pc_in = ch_callset.map { c, _bed, _bim, fam ->
        tuple(c, file(tpl(params.PhenoTpl, c), checkIfExists: true),
                 file(tpl(params.CovarTpl, c), checkIfExists: true), fam)
    }
    PREP_PHENO_COV(ch_pc_in, script_file('prep_pheno_cov.py'))
    ch_coding = PREP_PHENO_COV.out.pc.map { c, _rp, _rc, j -> tuple(c, j) }

    // ── STAGE 1 — the map ──────────────────────────────────────────────────
    SPLIT_BIM(ch_callset.map { c, _bed, bim, _fam -> tuple(c, bim) })
    ch_bim_chunks = byChrom(SPLIT_BIM.out.chunks, /bim\.chr([0-9]+)\.tsv/)

    MAKE_MINAC_KEEP(ch_callset.map { c, bed, bim, fam ->
        tuple(c, minac[c], bed, bim, fam) })

    MAP_CHUNK(ch_bim_chunks.combine(MAKE_MINAC_KEEP.out.keep, by: 0),
              script_file('build_gene_map.py'))

    // BARRIER B1 — the denominator is fixed here, before any P-value.
    ch_merge_map_in = MAP_CHUNK.out.map.groupTuple(by: 0)
        .join(MAP_CHUNK.out.anno.groupTuple(by: 0), by: 0)
        .join(MAP_CHUNK.out.qc.groupTuple(by: 0),   by: 0)
    MERGE_MAP(ch_merge_map_in, script_file('merge_gene_map.py'))

    // One tuple per (cohort, stratum), files matched by exact name.
    ch_stratum = MERGE_MAP.out.strata.flatMap { c, maps, idxs, vars, dens, cols ->
        params.Strata.collect { s ->
            tuple(c, s.tag,
                  pick(maps, 'map',            s.tag, '.tsv.gz'),
                  pick(idxs, 'gene_index',     s.tag, '.tsv'),
                  pick(vars, 'variants',       s.tag, '.txt'),
                  pick(dens, 'denominator',    s.tag, '.tsv'),
                  pick(cols, 'map_collisions', s.tag, '.tsv'))
        }
    }
    ch_gidx  = ch_stratum.map { c, s, _m, gi, _v, _d, _col -> tuple(c, s, gi) }
    ch_denom = ch_stratum.map { c, s, _m, _gi, _v, d, _col -> tuple(c, s, d) }
    ch_coll  = ch_stratum.map { c, s, _m, _gi, _v, _d, col -> tuple(c, s, col) }
    ch_map   = ch_stratum.map { c, s, m, _gi, _v, _d, _col -> tuple(c, s, m) }

    // ── STAGE 2 — test inputs ──────────────────────────────────────────────
    ch_emit_in = ch_stratum
        .map { c, s, m, gi, v, _d, _col -> tuple(c, s, m, gi, v) }
        .combine(ch_callset, by: 0)
        .combine(PREP_PHENO_COV.out.keep, by: 0)
        .combine(ch_coding, by: 0)
    EMIT_TEST_INPUTS(ch_emit_in, script_file('emit_test_inputs.py'))

    ch_setfile = EMIT_TEST_INPUTS.out.setfiles.flatMap { c, s, files ->
        (files instanceof List ? files : [files]).collect { f ->
            tuple(c, s, (f.name =~ /rvtest\.set\.chr([0-9]+)\.txt/)[0][1].toString(), f) } }

    // ── STAGE 3 — the associations ─────────────────────────────────────────
    ch_methods = channel.fromList(
        params.RvtestMethods.collect { mt -> [mt.tag, mt.opt, mt.out, (mt.wald ?: '')] })
    ch_pc_rvt  = PREP_PHENO_COV.out.pc.map { c, rp, rc, j -> tuple(c, rp, rc, j) }

    ch_rvtest_in = ch_setfile
        .combine(EMIT_TEST_INPUTS.out.vcf, by: [0, 1])
        .combine(ch_methods)
        .combine(ch_pc_rvt, by: 0)
        .map { c, s, chrom, setf, vcf, tbi, method, opt, sfx, wsfx, ph, cv, coding ->
               tuple(c, s, method, opt, sfx, wsfx, chrom, vcf, tbi, setf, ph, cv, coding) }
    RVTEST(ch_rvtest_in)

    // ── STAGE 4 — B2, B3, significance ─────────────────────────────────────
    ch_merge_in = RVTEST.out.assoc
        .map { c, s, m, chrom, a, w -> tuple(c, s, m, chrom, a, w) }
        .groupTuple(by: [0, 1, 2])
        .combine(ch_gidx, by: [0, 1])
        .combine(EMIT_TEST_INPUTS.out.setstats, by: [0, 1])
    MERGE_RVTEST(ch_merge_in, script_file('merge_rvtest.py'))

    ch_scan_in = MERGE_RVTEST.out.merged
        .map { c, s, _m, f -> tuple(c, s, f) }.groupTuple(by: [0, 1])
        .join(ch_gidx,  by: [0, 1])
        .join(ch_denom, by: [0, 1])
        .join(ch_coll,  by: [0, 1])
        .join(ch_map,   by: [0, 1])
    SCAN(ch_scan_in, script_file('scan.py'))

    ch_sig_in = SCAN.out.scan
        .map { c, s, _scan, qc -> tuple(c, s, qc) }
        .combine(ch_denom, by: [0, 1])
        .combine(channel.fromList(params.RvtestMethods.collect { mt -> mt.tag }))
        .map { c, s, qc, den, method -> tuple(c, s, method, qc, den) }
    SIGNIFICANCE(ch_sig_in, script_file('gene_significance.py'))

    // ── STAGE 5 — figures, per cohort ──────────────────────────────────────
    SCAN_SCALE(SCAN.out.scan.map { _c, _s, sc, _qc -> sc }.collect(),
               script_file('scan_scale.py'))
    // .first() makes the one scale file a VALUE channel: a queue channel holding
    // a single item would let only ONE of the six PLOT_SCAN tasks run.
    PLOT_SCAN(SCAN.out.scan.combine(ch_denom, by: [0, 1]),
              SCAN_SCALE.out.scale.first(),
              script_file('plot_gene_scan.py'))

    ch_anno_in = MERGE_MAP.out.anno
        .combine(ch_gidx.map  { c, _s, gi -> tuple(c, gi) }.groupTuple(by: 0), by: 0)
        .combine(ch_denom.map { c, _s, d  -> tuple(c, d)  }.groupTuple(by: 0), by: 0)
    PLOT_ANNOTATION(ch_anno_in, script_file('plot_annotation.py'))

    // ── STAGE 6 — across cohorts ───────────────────────────────────────────
    // BARRIER B4 — the 3 x 2 x 2 grid needs every cell at once.
    COLLECT_ALL(
        SCAN.out.scan.map { _c, _s, scan, _qc -> scan }.collect(),
        SCAN.out.scan.map { _c, _s, _scan, qc -> qc }.collect(),
        SIGNIFICANCE.out.sig.map { _c, _s, _m, f -> f }.collect(),
        SCAN.out.identity.map { _c, _s, f -> f }.collect(),
        ch_denom.map { _c, _s, d -> d }.collect(),
        EMIT_TEST_INPUTS.out.setstats.map { _c, _s, f -> f }.collect(),
        EMIT_TEST_INPUTS.out.varstats.map { _c, _s, f -> f }.collect(),
        EMIT_TEST_INPUTS.out.multi.map { _c, _s, f -> f }.collect(),
        ch_coding.map { _c, j -> j }.collect())

    // gene_index.<stratum>.tsv and anno_counts.tsv carry NO cohort or stratum
    // column -- their file name is the only key, and a stageAs rename would
    // destroy it. The keys are prepended HERE, where they still exist.
    ch_gidx_all = ch_gidx
        .map { c, s, gi ->
            def lines = gi.text.readLines().findAll { ln -> ln }
            def out = ["cohort\tstratum\t" + lines[0]]
            lines.drop(1).each { ln -> out << "${c}\t${s}\t${ln}".toString() }
            out.join('\n') + '\n'
        }
        .collectFile(name: 'gene_index_all.tsv', sort: true, keepHeader: true, skip: 1,
                     storeDir: "${params.out_dir}/_comparison/tables")
    MERGE_MAP.out.anno
        .map { c, anno, _qc ->
            def lines = anno.text.readLines().findAll { ln -> ln }
            def out = ["cohort\t" + lines[0]]
            lines.drop(1).each { ln -> out << "${c}\t${ln}".toString() }
            out.join('\n') + '\n'
        }
        .collectFile(name: 'annotation_all.tsv', sort: true, keepHeader: true, skip: 1,
                     storeDir: "${params.out_dir}/_comparison/tables")
        .subscribe { f -> log.info "assoc_gene | keyed annotation table -> ${f}" }

    ch_tables = COLLECT_ALL.out.tables
    ROBUSTNESS(ch_tables, ch_gidx_all, script_file('robust_genes.py'))

    PLOT_GRID(ch_tables.map { sc, qc, _sg, _id, _dn, _ss, _vs, _mc, _cs -> tuple(sc, qc) },
              script_file('plot_grid.py'))
    PLOT_ROBUST_GENES(
        ROBUSTNESS.out.robust.map { ev, ti, _vd, _md, cs -> tuple(ev, ti, cs) }
            .combine(ch_tables.map { _sc, _qc, _sg, _id, dn, _ss, _vs, _mc, _cs -> dn }),
        script_file('plot_robust_genes.py'))
    GENE_MODELS(ROBUSTNESS.out.windows,
                channel.fromPath("${params.shared_script_dir}/region_tracks.R"),
                channel.fromPath("${params.shared_script_dir}/gene_utils.R"))
    PLOT_GENE_DETAIL(
        ROBUSTNESS.out.robust.map { ev, ti, vd, md, cs -> tuple(ti, ev, vd, md, cs) },
        GENE_MODELS.out.models,
        script_file('plot_gene_detail.py'))

    // One README per fan-out family, from the members' stats.json files.
    ch_catalogue = PLOT_ANNOTATION.out.stats.collect()
        .map { fs -> tuple('variant_sets', '01.variant_sets', fs) }
        .mix(PLOT_SCAN.out.stats.collect().map { fs -> tuple('gene_scan', '02.gene_scan', fs) })
    CATALOGUE(ch_catalogue, script_file('catalogue.py'))

    // ── provenance ─────────────────────────────────────────────────────────
    def manifest = """{
  "component": "assoc_gene",
  "cohorts": ${groovy.json.JsonOutput.toJson(cohortList)},
  "cohort_nesting": "narrow_mainland is a strict subset of intermediate_mainland is a strict subset of full_mainland; the three differ almost entirely in controls (full_mainland has 20 more cases than narrow_mainland). Agreement across cohorts is a sample-selection sensitivity analysis, not replication.",
  "genotypes": "tuning.rv 02.callset_filter/filtered = 12_model_inputs fixed_model maf_lt_threshold after tuning.rv sample and variant QC",
  "minac": ${groovy.json.JsonOutput.toJson(minac)},
  "minac_source": "tuning.rv 05.qc_collect/minac_recommendation.tsv",
  "gene_membership": "snpEff gene field, ${params.GeneAssignment}; materialised as map.<stratum>.tsv.gz and handed to rvtest through --setFile as one 1-bp range per variant",
  "snpeff_index": "${params.SnpEffIndex}",
  "strata": ${groovy.json.JsonOutput.toJson(params.Strata.collect { st -> st.tag })},
  "strata_nesting": "moderate_high is a strict subset of low_moderate_high, by construction",
  "min_num_var": ${params.MinNumVar},
  "engine": "rvtest --setFile, fixed-effects logistic regression",
  "covariates": "${params.CovarName}",
  "methods": ${groovy.json.JsonOutput.toJson(params.RvtestMethods.collect { mt -> mt.tag })},
  "method_output_suffix": ${groovy.json.JsonOutput.toJson(params.RvtestMethods.collectEntries { mt -> [(mt.tag): mt.out] })},
  "effect_size": "OR and its 95% CI on cmc rows come from --burden cmcWald, a Wald fit of the SAME collapsed model whose score test produces the published p, run as a second model in the same rvtest invocation. z = beta/se does NOT reproduce the published p and the Wald p is deliberately not published. Kernel (skato) rows carry no effect size.",
  "skato_weights": "Beta(MAF; 1, 25), rvtest's default, asserted from the run log by RVTEST's post-flight gate; rho is the optimal SKAT-O mixing parameter rvtest reports",
  "missing_genotypes": "${params.ImputeMethod} (rvtest's default, passed explicitly and asserted from the log). Immaterial for CMC, whose collapse casts genotypes to int; used by the SKAT-O kernel.",
  "significance": "Per (cohort, stratum): Bonferroni alpha / n_genes_mapped with the denominator fixed at MERGE_MAP before any p-value exists (genes below min_num_var are excluded from the count), and Benjamini-Hochberg over exactly that family with genes rvtest did not return entering at p = 1. BH q < alpha is 'called'; Bonferroni is reported as the stricter tier. Because both rules share alpha and n, every Bonferroni call is a BH call (asserted in SCAN). No correction across methods or cohorts: they test the same genes on the same samples.",
  "alpha": ${params.Alpha},
  "robustness": "Tier 1 = called (BH) in every cohort and by both methods in at least one; Tier 2 = called in every cohort by at least one method; Tier 3 = called in at least two cohorts. Computed per (gene, stratum); Bonferroni never enters the rule. Descriptive, not a decision rule: no threshold, no p-value, no multiple-testing price.",
  "report_cohort": "${params.ReportCohort}",
  "carrier_mac_maf": "Per gene: carriers = samples with >= 1 minor-allele copy at any mapped variant, counted per group from a direct decode of tested.bed and cross-checked against plink2 --freq counts per variant; cumulative MAC per group from plink2 --freq counts; cumulative MAF = MAC / (2 * N_group). The minor allele is decided once on the pooled counts and applied to both groups, as rvtest does.",
  "span_overlap_rule": "rvtest's setFile carries one 1-bp range per variant, and a tabix interval matches a record over [POS, POS+len(REF)-1], so a multi-base-REF record is counted twice when a later mapped variant sits inside its span. scan.py predicts the over-count exactly and classifies such genes explained_by_span_overlap; the map identity stays an equality.",
  "strict_identity": ${params.StrictIdentity},
  "gene_models": "per-gene figures draw exon structure from Ensembl 86 (GRCh38) via the shared region_tracks.R; only the tested gene and genes overlapping it are drawn",
  "figures": "results/figures/0N.<family>/, five families in reading order, one README.md per family; figures draw lambda_GC only (the n_var split is tabulated); the main figure is 04.robust_genes/robust_genes.png"
}"""
    def index = """# assoc_gene results

| layer | where |
|---|---|
| deliverable: robust genes | `_comparison/robust_genes/summary_table.md` (+ `summary_table.<cohort>.tsv`, `tiers.tsv`, `evidence.tsv`, `variant_detail.tsv`, `multi_carrier_detail.tsv`) |
| main figure | `figures/04.robust_genes/robust_genes.png` |
| per-gene figures (Tier 1 and 2) | `figures/05.genes/<GENE>.png`, `index.tsv` |
| other figure families | `figures/01.variant_sets/`, `figures/02.gene_scan/`, `figures/03.calibration/` — one `README.md` each |
| cross-cohort tables | `_comparison/tables/` |
| per-cohort stages | `<cohort>/00.prep … 05.signals` |
| run record | `_run_info/run_manifest.json`, `trace.txt`, `report.html`, `timeline.html`, `dag.html` |

Cohorts: ${cohortList.join(', ')}. Variant sets: ${params.Strata.collect { st -> st.tag }.join(', ')}. Statistics: ${params.RvtestMethods.collect { mt -> mt.tag }.join(', ')}. Report cohort: ${params.ReportCohort}. Column dictionaries: `docs/OUTPUTS.md`; the figure system: `docs/FIGURES.md`."""
    WRITE_RUN_MANIFEST(channel.value(manifest), channel.value(index))
}
