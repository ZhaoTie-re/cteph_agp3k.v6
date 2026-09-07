#!/usr/bin/env nextflow
// =============================================================================
// hla.typing — HLA allele typing and amino-acid residue calling.
//
// CRAM -> HLA-region reads -> HLA-HD allele calls -> allele dosages and an
// IMGT-numbered residue matrix, for the full_mainland analysis cohort.
//
// THE SAMPLE SET IS THREE FILTERS, IN ORDER, AND THEY ARE NOT THE SAME KIND
//   Flag_JHRPv6 = True in the workbook          the v6 cohort      3,592
//   in wgs.auto.par's sample_qc.keep.id         analysable         3,569  (-23)
//   in PopGMM_output/full_mainland.fid_iid.txt  COMPARABLE         3,101  (-468)
//
//   The third is an ANCESTRY filter, not a quality one, and it is the reason this
//   component's scope is a cohort rather than "everything that typed". The only
//   external check available here asks whether an allele we called appears in a
//   mainland Japanese reference panel; a control who is not mainland Japanese
//   fails that check for a reason that has nothing to do with typing. Typing the
//   whole QC-passing set and subsetting downstream is the other valid design. This
//   one is chosen so every table published here has ONE denominator.
//
//   The 468 samples the cohort filter removes were typed by an earlier run and are
//   NOT deleted: their reads, HLA-HD output and the 3,569-sample tables built from
//   them are kept under results/_superseded.3569/, which is regenerable from work/.
//
// WHAT THIS COMPONENT DOES NOT DECIDE
//   Nothing about association. It produces genotypes; whether a residue is
//   associated with anything is a question for the component that reads them.
//   It also does not decide the sample set: the three files above do.
//
// THE THREE THINGS THAT ARE NOT NEGOTIABLE HERE
//
//   1. HLA-HD cannot report failure. Its estimation step is a generated shell
//      script with no `set -e` that checks none of its children, so an
//      out-of-memory kill of hla_est — the common failure, and it takes HLA-A and
//      HLA-B first — still exits 0 and publishes a truncated result. The truncated
//      line carries NO gene name, so a parser reading by row position silently
//      re-labels every gene after it.
//
//      validate_typing.py is the only signal there is, and it runs INSIDE the HLAHD
//      process, as the last line of its script. That placement is the whole point:
//      a truncation becomes a non-zero exit of the typing task itself, which is
//      what the memory ladder retries. It used to be a separate downstream process,
//      where its failure could only terminate the run — a downstream process cannot
//      re-run an upstream one — and the truncated result had already been published
//      by the time it looked at it.
//
//   2. The index is not beside the CRAM. 3,117 of 3,569 CRAM paths are symlinks
//      whose .crai sits beside the real file. The index is resolved by realpath in
//      BUILD_MANIFEST and passed explicitly to `samtools view -X`.
//
//   3. Residues come from IMGT's published alignments, never from a recomputed
//      one. Re-deriving them (translate exons, blastx, MUSCLE) replaces IMGT
//      numbering with an alignment column index that restarts per exon and moves
//      with every release, so a result at "DRB1 position 13" is neither
//      reproducible nor the position anyone else means.
//
// SIZE, AND WHY EVERYTHING BIG IS A SYMLINK
//   Per sample: ~90 MB of extracted FASTQ and ~1.2 GB of HLA-HD output, of which
//   98 % is intermediate SAM. Across the cohort that is 321 GB + 4.2 TB. All of it
//   is kept and NONE of it is copied — the tables, figures and run record that a
//   reader actually opens are copied, everything else is a symlink into work/.
//   The consequence is explicit: DO NOT CLEAN work/. Doing so leaves the tables
//   intact and removes every intermediate, and re-deriving them means re-reading
//   88 TB of CRAM.
//
// Layout
//   results/00.manifest/   sample_manifest.tsv · contig_list.txt · crai/
//           01.reads/      <sample>.hla.R{1,2}.fastq.gz          (symlinked)
//           02.typing/     <sample>/{result,log,mapfile,exon,…}  (symlinked)
//           03.alleles/    allele_calls · allele_ambiguity · typing_status
//           04.residues/   reference (symlinked) · allele_dosage · residue
//                          diplotype/dosage/sites · allele_match_map · unmatched
//           05.qc/         typing_qc_{sample,platform,group,gene}.tsv ·
//                          allele_frequency_{check,summary,sample}[.1kg].tsv ·
//                          allele_pgroup_map[.1kg].tsv · typing_confound.tsv
//           figures/       typing_qc · allele_frequency · allele_frequency_loci ·
//                          typing_confound — each .png with its .md
//           _run_info/     trace · report · timeline · dag · run_manifest · facts.json
//           _superseded.3569/  the 468 samples outside the cohort, and the tables
//                          built before the cohort filter existed. Not an output.
// =============================================================================

nextflow.enable.dsl = 2

// =============================================================================
// SITE CONFIGURATION — everything below is environment- or project-specific.
// Override on the command line (`--key value`) or with `-params-file cfg.yaml`;
// nothing outside this block needs to change to run the component elsewhere.
// =============================================================================
params.project_dir       = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6'
params.base_dir          = "${params.project_dir}/hla.typing"
params.script_dir        = "${params.base_dir}/scripts"
params.shared_script_dir = "${params.project_dir}/analysis/_shared/scripts"
params.out_dir           = "${params.base_dir}/results"

// ---- Sample sheet ----
params.sample_info       = "${params.project_dir}/info/cteph_agp3k.v6.20260507.xlsx"
params.IdCol             = 'ID_JHRPv6'
params.FlagCol           = 'Flag_JHRPv6'
params.CramCol           = 'Cram_Path'
params.PlatformCol       = 'WGS_Platform'
params.GroupCol          = 'Outcome'
params.SexCol            = 'Sex'
params.DepthCol          = 'Observed_Depth'
// The workbook flag defines the v6 cohort; wgs.auto.par's sample QC decides which
// of it is analysable. This component inherits that decision instead of making its
// own — the 23 exclusions are heterozygosity outliers, and a sample whose
// heterozygosity is anomalous is exactly the one least likely to resolve into two
// clean HLA haplotypes.
params.sample_qc_keep    = "${params.project_dir}/wgs.auto.par/results/07_sample_qc/run_qc/cteph_agp3k_v6_wgs_merged.sample_qc.keep.id"
// The analysis cohort. Applied after sample QC and kept separate from it: the two
// answer different questions and collapsing them would lose the distinction. It is
// a fixed list on disk — the union of 17 components of a Gaussian mixture fitted in
// BBJ principal-component space; the clustering itself is not this repository's code
// and PopGMM_output/keep_list_summary.tsv is its only in-repo provenance. Set to ''
// to type every QC-passing sample instead.
params.cohort_keep       = "${params.project_dir}/PopGMM_output/full_mainland.fid_iid.txt"
params.CohortName        = 'full_mainland'
params.MaxSamples        = 0        // 0 = all; a positive value takes a pilot
                                    // spread across every platform

// ---- Reference resources ----
// Our own hs38DH, the one the CRAMs were aligned to: 3,366 sequences, matching
// the CRAM @SQ exactly, including the 525 HLA* decoy contigs typing depends on.
params.reference         = "${params.project_dir}/../nagasaki_pipeline/data/hs38DH.fa"
params.MhcChrom          = 'chr6'
params.MhcStart          = 29540000
params.MhcEnd            = 33420000

// ---- HLA-HD ----
// The dictionary and the IMGT alignments are PINNED COPIES. The directories they
// came from are owned by other users and change without notice; a reference that
// moves under a study is not a reference. See software/hlahd.dictionary/README.md.
params.hlahd_dir         = '/LARGE0/gr10478/bin/hlahd.1.7.1'
params.hlahd_sh          = "${params.hlahd_dir}/bin/hlahd.sh"
params.hlahd_freq        = "${params.hlahd_dir}/freq_data"
// Pinned beside the dictionary, NOT read from the group install: this 1 KB file
// defines the gene list every result is validated against, so a silent swap of it
// would start failing samples whose typing was fine.
params.hlahd_split       = '/LARGE0/gr10478/b37974/software/hlahd.dictionary/HLA_gene.split.3.50.0.txt'
params.hlahd_dict        = '/LARGE0/gr10478/b37974/software/hlahd.dictionary/IMGT-HLA-3.64.0_N150'
params.imgt_alignments   = '/LARGE0/gr10478/b37974/software/hlahd.dictionary/IMGT-alignments-3.64.0'
params.MinReadLength     = 100      // hlahd.sh -m
params.CuttingRate       = 0.95     // hlahd.sh -c

// ---- Genes ----
// HLA-HD types 33. IMGT publishes a protein alignment for 19 of them, and only
// those can carry residues. The other 14 keep their allele calls and are recorded
// as having no residue reference rather than being dropped without a word.
params.ResidueGenes      = 'A,B,C,DMA,DMB,DOA,DOB,DPA1,DPB1,DQA1,DQB1,DRA,DRB1,DRB3,DRB4,DRB5,E,F,G'
params.MinResidueCount   = 1        // drop a residue (or allele) seen on fewer
                                    // chromosomes than this
params.AlleleFieldDepth  = 2        // the depth allele dosage is encoded at, and
                                    // the depth the external reference carries

// ---- External frequency reference ----
// The one external check this component can afford. See jMorp_HLA_types/README.md and
// docs/OPEN_QUESTIONS.md §2 for what it does and does not establish.
//
// jMorp 61KJPN-HLA is 61,424 Japanese individuals over 13 loci and is PRIMARY. The
// 1000 Genomes JPT panel is 105 individuals over 5 loci, missing DPB1, DQA1 and
// DRB3/4/5 -- the loci this component types least reliably -- and is SECONDARY. It is
// not a fallback: both run on every execution, and the pair is what licenses reading
// the group contrast rather than the level. An allele absent from 105 people is weak
// evidence; absent from 61,424 Japanese it is strong, which is why the panel a caveat
// names has to be the panel that produced the figure.
params.hla_truth_jmorp   = "${params.project_dir}/../jMorp_HLA_types/HLA_allele_frequencies_61K.txt"
params.hla_truth_1kg     = "${params.project_dir}/../1KG_HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt"
// BOTH panels are run, every time. jMorp is primary and 1000G JPT is the secondary
// check, and running both is not redundancy: on identical calls the two disagree by
// more than a factor of two on the ABSOLUTE unconfirmed rate while agreeing closely
// on the RATIO between our two groups. That is the demonstration — not the assertion
// — that the level is a property of the panel and only the contrast is readable.

// IPD-IMGT/HLA P-group definitions, pinned to the SAME 3.64.0 release as
// params.imgt_alignments so the alignments and the groupings cannot drift apart.
// jMorp names alleles by P group and HLA-HD does not; without this, DRB4*01:03 -- 76 %
// of our DRB4 chromosomes -- reads as absent from a 61,424-person panel. Used ONLY by
// ALLELE_FREQ_CHECK. It never reaches allele_dosage.tsv or residue_dosage.tsv, which
// the association reads: P groups merge alleles that differ outside the antigen
// recognition domain, and the residue analysis reads the full protein.
params.hla_pgroups       = "${params.project_dir}/../jMorp_HLA_types/hla_nom_p.txt"

// How the reference is named on the figure, and which loci carry a caveat that has to be
// visible ON the panel rather than only in the document beside it.
//
// FlagLoci is a LABEL, not the decision. Whether a locus is comparable at all is decided
// in allele_freq_check.py from the two denominators, mechanically, and written to the
// `comparable` column: jMorp assigns a DRB3 allele to every chromosome in the panel while
// this pipeline calls DRB3 on about half of them, so the two are not measuring the same
// quantity. See OPEN_QUESTIONS.md section 3.
params.TruthName         = 'jMorp 61KJPN-HLA'
params.TruthCite         = 'ToMMo jMorp; Tadaka et al. 2023'
params.TruthName1kg      = '1000 Genomes JPT'
params.TruthCite1kg      = '1000 Genomes Project; Abi-Rached et al. 2018'
params.FlagLoci          = 'DRB3'

params.TruthPopulation   = 'JPT'    // 1kg_wide only; the only Japanese population there
params.ControlGroup      = 'AGP3K'  // the `group` value that is a population sample
// ONE LOCUS LIST PER PANEL, and they must not be shared. A locus the panel does not
// carry scores as entirely unconfirmed rather than as absent, so running the 13-locus
// list against the 5-locus panel reports 62.8 % unconfirmed instead of 10.4 % — a
// number about the list, not about the typing. Measured, not hypothetical.
params.RefLociJmorp      = 'A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1,DRB3,DRB4,E,F,G'  // the 13 jMorp carries
params.RefLoci1kg        = 'A,B,C,DQB1,DRB1'                                 // the 5 the 1KG panel has

// ---- Environment ----
params.conda_env         = 'cteph_geno_pro'
params.samtools          = '/LARGE0/gr10478/b37974/software/samtools/bin/samtools'
params.bowtie2_dir       = '/LARGE0/gr10478/b37974/software/bowtie2-2.5.5-linux-x86_64'

// =============================================================================
// HELPERS
// =============================================================================

// Scripts are staged as path inputs, never referenced as bare param paths: only a
// path input is content-hashed, so this is what makes editing a script actually
// invalidate its task on -resume instead of silently reusing the old result.
def script_file(String name) {
    return file("${params.script_dir}/${name}", checkIfExists: true)
}

// =============================================================================
// PROCESSES
// =============================================================================

/* STEP 1 · BUILD_MANIFEST  -> 00.manifest/ — the sample set, with every CRAM and index resolved
 * and verified before a single job is submitted. Aborts naming the samples.
 *
 * IT DOES NOT LOOK IN 00.manifest/crai/, AND THAT IS DELIBERATE. build_manifest.py can
 * (--crai-dir), and wiring it here was tried: it lets a run whose work/ was lost reuse the
 * 111 indexes INDEX_CRAM already built, saving about 24 CPU-hours. The cost is worse than
 * the saving — the manifest then depends on this pipeline's OWN previous output, so it is
 * not idempotent, and the second run writes a different `crai` column for those 111 samples,
 * which re-runs their extraction and their typing. Measured: a -resume that should have been
 * free instead re-ran 111 of 3,569 EXTRACT_READS and HLAHD tasks. The manifest describes the
 * INPUTS, and it has to be a function of the inputs alone. */
process BUILD_MANIFEST {
    executor 'local'
    tag 'manifest'
    publishDir "${params.out_dir}/00.manifest", mode: 'copy'

    input:
    path script
    path xlsx
    path keep_id
    path cohort_keep

    output:
    path('sample_manifest.tsv'), emit: manifest

    script:
    // An empty cohort list is staged as the sentinel NO_COHORT, which does not
    // exist as a file, so the flag is dropped rather than pointed at nothing.
    def cohort_arg = cohort_keep.name == 'NO_COHORT' ? '' : "--cohort-keep ${cohort_keep}"
    """
    source activate ${params.conda_env}
    python3 ${script} \\
        --xlsx ${xlsx} \\
        --id-col ${params.IdCol} --flag-col ${params.FlagCol} \\
        --cram-col ${params.CramCol} --platform-col ${params.PlatformCol} \\
        --group-col ${params.GroupCol} --sex-col ${params.SexCol} \\
        --depth-col ${params.DepthCol} \\
        --keep-id ${keep_id} \\
        ${cohort_arg} \\
        --max-samples ${params.MaxSamples} \\
        --out sample_manifest.tsv
    """
}

/* STEP 2 · CONTIG_LIST  -> 00.manifest/ — the extraction regions, built once for the run from the
 * .fai. Reading them from the 3 GB FASTA in every task, as the obvious
 * implementation does, costs ~200 s per sample and buys nothing. */
process CONTIG_LIST {
    executor 'local'
    tag 'contigs'
    publishDir "${params.out_dir}/00.manifest", mode: 'copy'

    input:
    path script
    path fai

    output:
    path('contig_list.txt'), emit: contigs

    script:
    """
    source activate ${params.conda_env}
    python3 ${script} --fai ${fai} \\
        --mhc-chrom ${params.MhcChrom} \\
        --mhc-start ${params.MhcStart} --mhc-end ${params.MhcEnd} \\
        --out contig_list.txt
    """
}

/* STEP 3 · INDEX_CRAM  -> 00.manifest/crai/ — build the index for a CRAM that has none, into our own
 * tree. 111 of this study's CRAMs are unindexed; nothing is ever written next to
 * the source, which lives in a directory we do not own. */
process INDEX_CRAM {
    executor 'slurm'
    queue    'gr10478b'
    time     '6h'
    tag      "${sample}"
    publishDir "${params.out_dir}/00.manifest/crai", mode: 'copy'

    input:
    tuple val(sample), path(cram)

    output:
    tuple val(sample), path("${sample}.cram.crai"), emit: crai

    script:
    """
    set -euo pipefail
    ${params.samtools} index -@ ${task.ext.threads ?: 1} -o ${sample}.cram.crai ${cram}
    """
}

/* STEP 4 · EXTRACT_READS  -> 01.reads/ — CRAM to the FASTQ pair HLA-HD reads. */
process EXTRACT_READS {
    executor 'slurm'
    queue    'gr10478b'
    time     '4h'
    tag      "${sample}"

    // Symlinked, not copied. Producing these is the most expensive step in the
    // component — it reads 88 TB of CRAM for ~360 CPU-hours — and they are only
    // ~90 MB per sample. Anything that re-types (a newer IMGT dictionary, a
    // different -m, a second tool for validation) starts from here rather than
    // from the CRAMs. They are worth being able to FIND, which a work-dir hash
    // does not allow.
    publishDir "${params.out_dir}/01.reads", mode: 'symlink'

    input:
    tuple val(sample), path(cram), path(crai)
    path script
    path contigs
    path reference
    path ref_fai

    output:
    tuple val(sample), path("${sample}.hla.R1.fastq.gz"),
                       path("${sample}.hla.R2.fastq.gz"), emit: fastq
    path("${sample}.hla.singleton.fastq.gz"),             emit: singleton

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} \\
        --sample ${sample} \\
        --cram ${cram} --crai ${crai} \\
        --reference ${reference} \\
        --contig-list ${contigs} \\
        --samtools ${params.samtools} \\
        --threads ${task.ext.threads ?: 1}
    """
}

/* STEP 5 · HLAHD  -> 02.typing/ — typing, AND the gate on it. The whole output tree
 * is kept, including the alignments: they are 98 % of 1.2 GB per sample, which is why
 * it is symlinked rather than copied, but they are the only record of how a call was
 * reached. */
process HLAHD {
    executor 'slurm'
    queue    'gr10478b'
    time     '8h'
    tag      "${sample}"
    // Symlinked, and the alignment tree is KEPT. It is 1.17 GB per sample —
    // 4.2 TB across the cohort — which is why copying it is out of the question,
    // but a symlink costs nothing and the SAMs are the only record of how each
    // allele call was reached. Re-deriving them means re-running bowtie2 over the
    // whole cohort.
    //
    // publishDir runs only for a task that SUCCEEDED, and success here now includes
    // passing validation. Nothing truncated reaches 02.typing.
    publishDir "${params.out_dir}/02.typing", mode: 'symlink'

    input:
    tuple val(sample), path(fq1), path(fq2)
    path validator
    path gene_split

    output:
    tuple val(sample), path("${sample}"), emit: typed

    script:
    """
    set -euo pipefail
    export PATH=${params.hlahd_dir}/bin:${params.bowtie2_dir}:\$PATH
    mkdir -p ${sample}

    # hlahd.sh returns 0 whatever happens inside it; the exit code below is
    # therefore not evidence of anything. The validator at the end of this script is.
    ${params.hlahd_sh} \\
        -t ${task.ext.threads ?: 1} \\
        -m ${params.MinReadLength} \\
        -c ${params.CuttingRate} \\
        -f ${params.hlahd_freq} \\
        ${fq1} ${fq2} \\
        ${params.hlahd_split} \\
        ${params.hlahd_dict} \\
        ${sample} ./ \\
        > hlahd.run.log 2>&1 || true

    # The uncompressed FASTQ copies HLA-HD leaves behind are the one thing not
    # worth keeping: they are a duplicate of 01.reads, uncompressed.
    find ${sample} -name '*.fastq' -delete || true
    mkdir -p ${sample}/result ${sample}/log
    # Inside the sample directory, beside the per-gene logs it belongs with —
    # not a sibling of it, which put two entries per sample in one flat directory.
    mv hlahd.run.log ${sample}/hlahd.run.log

    # THE LAST LINE, and it must stay last. A truncated result exits non-zero HERE,
    # in the task the memory ladder retries, rather than in a downstream process that
    # can only kill the run. The per-sample record is written beside the result.
    source activate ${params.conda_env}
    python3 ${validator} \\
        --sample ${sample} \\
        --result-dir ${sample}/result \\
        --run-log ${sample}/hlahd.run.log \\
        --gene-split ${gene_split} \\
        --out ${sample}/typing_validation.tsv
    """
}

/* STEP 6 · COLLECT_ALLELES  -> 03.alleles/ — every result into one table, parsed by gene name. */
process COLLECT_ALLELES {
    executor 'slurm'
    queue    'gr10478b'
    time     '2h'
    tag      'collect'
    publishDir "${params.out_dir}/03.alleles", mode: 'copy'

    // Staged under HLA-HD's own names, which are already unique per sample. The
    // earlier `stageAs: 'res_*'` needed a parallel list of `SAMPLE=res_N` arguments
    // built in the workflow, so two independent orderings had to agree — and
    // Nextflow does not number `res_*` at all when the collection holds one file,
    // which made a single-sample run fail outright.
    input:
    path results
    val n_expected
    path script
    path gene_split

    output:
    path('allele_calls.tsv'),     emit: calls
    path('allele_ambiguity.tsv'), emit: ambiguity
    path('typing_status.tsv'),    emit: status

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} \\
        --result-dir . \\
        --expect-n ${n_expected} \\
        --gene-split ${gene_split} \\
        --out-calls allele_calls.tsv \\
        --out-ambiguity allele_ambiguity.tsv \\
        --out-status typing_status.tsv
    """
}

/* STEP 7 · BUILD_RESIDUE_REF  -> 04.residues/ — IMGT's alignments to a residue reference. */
process BUILD_RESIDUE_REF {
    executor 'slurm'
    queue    'gr10478b'
    time     '2h'
    tag      'residue-ref'
    // The reference is 34 MB and is identical for every run against the same IMGT
    // release — it is derived data, not a result of this cohort. Symlinked; the
    // summary that describes it is copied, so the record survives work/ being cleaned.
    publishDir "${params.out_dir}/04.residues", mode: 'symlink', pattern: 'residue_reference'
    publishDir "${params.out_dir}/04.residues", mode: 'copy',    pattern: '*.tsv'

    input:
    path script
    path alignments

    output:
    path('residue_reference'),              emit: ref
    path('residue_reference_summary.tsv'),  emit: summary

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} \\
        --alignment-dir ${alignments} \\
        --genes ${params.ResidueGenes} \\
        --out-dir residue_reference \\
        --out-summary residue_reference_summary.tsv
    """
}

/* STEP 8 · RESIDUE_MATRIX  -> 04.residues/ — calls x reference to allele dosages,
 * residue diplotypes and residue dosages. */
process RESIDUE_MATRIX {
    executor 'slurm'
    queue    'gr10478b'
    time     '4h'
    tag      'residues'
    publishDir "${params.out_dir}/04.residues", mode: 'copy'

    input:
    path calls
    path reference
    path script

    output:
    path('allele_dosage.tsv'),     emit: allele_dosage
    path('residue_diplotype.tsv'), emit: diplotype
    path('residue_dosage.tsv'),    emit: dosage
    path('residue_sites.tsv'),     emit: sites
    path('allele_match_map.tsv'),  emit: match_map
    path('allele_unmatched.tsv'),  emit: unmatched

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} \\
        --calls ${calls} --reference-dir ${reference} \\
        --genes ${params.ResidueGenes} \\
        --min-count ${params.MinResidueCount} \\
        --allele-field-depth ${params.AlleleFieldDepth} \\
        --out-allele-dosage allele_dosage.tsv \\
        --out-diplotype residue_diplotype.tsv \\
        --out-dosage residue_dosage.tsv \\
        --out-sites residue_sites.tsv \\
        --out-match-map allele_match_map.tsv \\
        --out-unmatched allele_unmatched.tsv
    """
}

/* STEP 9 · TYPING_QC  -> 05.qc/ — quality per platform, the axis that matters here,
 * and per case/control group, which a reader asks for and which is confounded. */
process TYPING_QC {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      'qc'
    publishDir "${params.out_dir}/05.qc", mode: 'copy', pattern: '*.tsv'

    input:
    tuple path(calls), path(ambiguity), path(manifest)
    path script

    output:
    path('typing_qc_sample.tsv'),   emit: sample_qc
    path('typing_qc_platform.tsv'), emit: platform_qc
    path('typing_qc_group.tsv'),    emit: group_qc
    path('typing_qc_gene.tsv'),     emit: gene_qc

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} \\
        --calls ${calls} --ambiguity ${ambiguity} \\
        --manifest ${manifest} --genes ${params.ResidueGenes} \\
        --out-sample typing_qc_sample.tsv \\
        --out-platform typing_qc_platform.tsv \\
        --out-group typing_qc_group.tsv \\
        --out-gene typing_qc_gene.tsv
    """
}

/* STEP 10 · ALLELE_FREQ_CHECK  -> 05.qc/ — the only external check here: do the
 * control allele frequencies reproduce the published Japanese ones? Everything else
 * in 05.qc measures whether a call was MADE, and none of that can see a call that was
 * made confidently and is wrong.
 *
 * Also emits allele_pgroup_map.tsv, the crosswalk the association reads: every allele
 * this cohort carries, its P group, and the reference frequency behind it. It is
 * emitted HERE and not from RESIDUE_MATRIX on purpose — touching that process would
 * regenerate allele_dosage.tsv and residue_dosage.tsv, which the association reads and
 * which must not move because a QC reference changed. */
process ALLELE_FREQ_CHECK {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      'freq-check'
    publishDir "${params.out_dir}/05.qc", mode: 'copy'

    input:
    tuple val(meta), path(truth), path(calls), path(manifest)
    path script
    path pgroups

    output:
    tuple val(meta), path("allele_frequency_check${meta.suffix}.tsv"),
                     path("allele_frequency_summary${meta.suffix}.tsv"),
                     path("allele_frequency_sample${meta.suffix}.tsv"),
                     path("allele_pgroup_map${meta.suffix}.tsv"), emit: out

    script:
    // --population is meaningful only for the wide 1KG layout, which has one column
    // per population; the jMorp file is already one population and rejects it.
    def pop_arg = meta.population ? "--population ${meta.population}" : ''
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} \\
        --calls ${calls} --manifest ${manifest} \\
        --truth ${truth} \\
        --truth-format ${meta.format} \\
        --p-groups ${pgroups} \\
        ${pop_arg} \\
        --control-group ${params.ControlGroup} \\
        --genes ${meta.loci} \\
        --field-depth ${params.AlleleFieldDepth} \\
        --out allele_frequency_check${meta.suffix}.tsv \\
        --out-summary allele_frequency_summary${meta.suffix}.tsv \\
        --out-sample allele_frequency_sample${meta.suffix}.tsv \\
        --out-pgroup-map allele_pgroup_map${meta.suffix}.tsv
    """
}

/* STEP 11 · PLOT_TYPING_QC  -> figures/ — one standalone figure, so it carries its
 * own document rather than a family catalogue. */
process PLOT_TYPING_QC {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      'qc-figures'
    publishDir "${params.out_dir}/figures", mode: 'copy'

    input:
    tuple path(sample_qc), path(gene_qc), path(sites)
    path script

    output:
    path('*.png'), emit: png
    path('*.md'),  emit: doc

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} \\
        --sample-qc ${sample_qc} \\
        --gene-qc ${gene_qc} --sites ${sites} \\
        --out-png typing_qc.png
    """
}

/* STEP 12 · PLOT_ALLELE_FREQ  -> figures/ — the frequency check, drawn, as TWO
 * figures from one script. `allele_frequency.png` answers "do the frequencies
 * agree?" and `allele_frequency_loci.png` is the per-locus record behind it. They
 * were one 17.4-inch figure in which 13 near-identical scatters crowded out the
 * summary they were supposed to support; a reader who wants the answer and a reader
 * who wants the evidence are not served by the same panel. */
process PLOT_ALLELE_FREQ {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      'freq-figure'
    publishDir "${params.out_dir}/figures", mode: 'copy'

    input:
    tuple path(check), path(summary), path(summary_2nd)
    path script

    output:
    path('*.png'), emit: png
    path('*.md'),  emit: doc

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} \\
        --check ${check} --summary ${summary} \\
        --summary-secondary ${summary_2nd} \\
        --reference-name '${params.TruthName}' \\
        --reference-name-secondary '${params.TruthName1kg}' \\
        --flag-genes '${params.FlagLoci}' \\
        --out-png allele_frequency.png \\
        --out-png-loci allele_frequency_loci.png
    """
}

/* STEP 13 · PLOT_TYPING_CONFOUND  -> figures/ — the one quality problem this
 * component has, drawn. typing_qc.png measures completeness, which is saturated
 * here; this measures whether the two groups are typed EQUALLY WELL, and answers the
 * obvious follow-up (is it depth?) with the data.
 *
 * It reads BOTH panels' per-sample tables. The absolute unconfirmed rate is a
 * property of the panel and the contrast between groups is not, and a figure that
 * shows one panel can only assert that in its caption. Showing both makes the
 * argument visible: the level moves, the ratio does not. */
process PLOT_TYPING_CONFOUND {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      'confound-figure'
    publishDir "${params.out_dir}/figures", mode: 'copy', pattern: '*.png'
    publishDir "${params.out_dir}/figures", mode: 'copy', pattern: '*.md'
    publishDir "${params.out_dir}/05.qc",   mode: 'copy', pattern: '*.tsv'

    input:
    tuple path(sample), path(sample_2nd)
    path script

    output:
    path('typing_confound.png'), emit: png
    path('typing_confound.md')
    path('typing_confound.tsv')

    script:
    """
    set -euo pipefail
    source activate ${params.conda_env}
    python3 ${script} \\
        --sample ${sample} --reference-name '${params.TruthName}' \\
        --sample-secondary ${sample_2nd} --reference-name-secondary '${params.TruthName1kg}' \\
        --control-group ${params.ControlGroup} \\
        --out-png typing_confound.png --out-table typing_confound.tsv
    """
}

/* STEP 14 · WRITE_RUN_MANIFEST  -> _run_info/ — what actually ran. */
process WRITE_RUN_MANIFEST {
    executor 'local'
    tag 'manifest'
    publishDir "${params.run_info_dir}", mode: 'copy'

    input:
    val manifest_json

    output:
    path('run_manifest.json')

    script:
    """
    cat > run_manifest.json <<'MANIFEST'
${manifest_json}
MANIFEST
    """
}

// =============================================================================
// WORKFLOW
// =============================================================================
workflow {

    // -- [1,2,3] the sample set and its regions, verified before anything runs ---
    // An absent cohort list is a real configuration — type every QC-passing sample —
    // so it is a sentinel name rather than a missing file, which `path` cannot stage.
    ch_cohort = params.cohort_keep
        ? file(params.cohort_keep, checkIfExists: true)
        : file("${projectDir}/NO_COHORT")

    BUILD_MANIFEST(script_file('build_manifest.py'),
                   file(params.sample_info, checkIfExists: true),
                   file(params.sample_qc_keep, checkIfExists: true),
                   ch_cohort)

    CONTIG_LIST(script_file('contig_list.py'),
                file("${params.reference}.fai", checkIfExists: true))

    ch_sample = BUILD_MANIFEST.out.manifest
        .splitCsv(header: true, sep: '\t')
        .map { r -> tuple(r.sample_id, file(r.cram), r.crai) }

    // A blank crai means no index exists anywhere; build one into our own tree.
    ch_have  = ch_sample.filter { _s, _c, i -> i     }.map { s, c, i -> tuple(s, c, file(i)) }
    ch_need  = ch_sample.filter { _s, _c, i -> !i    }.map { s, c, _i -> tuple(s, c) }
    ch_built = INDEX_CRAM(ch_need).crai
        .join(ch_need.map { s, c -> tuple(s, c) })
        .map { s, i, c -> tuple(s, c, i) }

    ch_ready = ch_have.mix(ch_built)

    // -- [4] extract, [5] type and validate in one task -------------------------
    // The contig list is consumed once per sample, which works because
    // `CONTIG_LIST.out.contigs` is a VALUE channel: a process emits one when all of
    // its inputs are values, and both of CONTIG_LIST's are `file()`. A value channel
    // can be read an unlimited number of times. Give CONTIG_LIST a queue input and
    // that stops holding — extraction would then run for ONE sample, not 3,569 —
    // so keep its inputs values.
    EXTRACT_READS(ch_ready,
                  script_file('extract_hla_reads.py'),
                  CONTIG_LIST.out.contigs,
                  file(params.reference, checkIfExists: true),
                  file("${params.reference}.fai", checkIfExists: true))

    // The validator runs inside HLAHD, so `typed` carries only samples that passed
    // it. A sample that did not was retried with more memory and, if it still
    // failed, terminated the run — it is never silently absent from what follows.
    HLAHD(EXTRACT_READS.out.fastq,
          script_file('validate_typing.py'),
          file(params.hlahd_split, checkIfExists: true))

    // -- [6] one allele table ---------------------------------------------------
    // Only `<sample>_final.result.txt` is passed on, never the directory it sits in:
    // that directory holds 67 files per sample, and staging all of them into the
    // single COLLECT_ALLELES task would be 239,123 symlinks for this cohort against
    // the 3,569 it actually needs.
    ch_result = HLAHD.out.typed
        .map { s, d -> file("${d}/result/${s}_final.result.txt") }
        .toSortedList { a, b -> a.name <=> b.name }

    COLLECT_ALLELES(ch_result,
                    ch_result.map { rows -> rows.size() },
                    script_file('collect_alleles.py'),
                    file(params.hlahd_split, checkIfExists: true))

    // -- [7,8] residues and allele dosages --------------------------------------
    BUILD_RESIDUE_REF(script_file('build_residue_ref.py'),
                      file(params.imgt_alignments, checkIfExists: true))

    RESIDUE_MATRIX(COLLECT_ALLELES.out.calls,
                   BUILD_RESIDUE_REF.out.ref,
                   script_file('residue_matrix.py'))

    // -- [9,10] QC, internal then external --------------------------------------
    // typing_status.tsv is not staged here: typing_qc.py derives every outcome
    // count from allele_calls.tsv itself, and the flag that used to pass it was
    // declared required and never read. verify.sh section 4 is this rule.
    TYPING_QC(COLLECT_ALLELES.out.calls
                  .combine(COLLECT_ALLELES.out.ambiguity)
                  .combine(BUILD_MANIFEST.out.manifest),
              script_file('typing_qc.py'))

    // Both panels, one task each. The per-panel locus list travels WITH the panel:
    // sharing one list between them is what silently reported 62.8 % unconfirmed
    // against the 5-locus 1KG file, because the 8 loci it does not carry counted as
    // unconfirmed rather than as absent.
    ch_panel = channel.of(
        [[tag: 'jmorp', role: 'primary',   format: 'jmorp_long', suffix: '',
          loci: params.RefLociJmorp, population: '',
          name: params.TruthName,    cite: params.TruthCite],
         file(params.hla_truth_jmorp, checkIfExists: true)],
        [[tag: '1kg',   role: 'secondary', format: '1kg_wide',   suffix: '.1kg',
          loci: params.RefLoci1kg,   population: params.TruthPopulation,
          name: params.TruthName1kg, cite: params.TruthCite1kg],
         file(params.hla_truth_1kg,   checkIfExists: true)])

    ALLELE_FREQ_CHECK(ch_panel.combine(COLLECT_ALLELES.out.calls
                                           .combine(BUILD_MANIFEST.out.manifest)),
                      script_file('allele_freq_check.py'),
                      file(params.hla_pgroups, checkIfExists: true))

    ch_freq = ALLELE_FREQ_CHECK.out.out.branch { meta, _c, _s, _sm, _pg ->
        primary:   meta.role == 'primary'
        secondary: true
    }

    // -- [11,12] figures --------------------------------------------------------
    // typing_qc_platform.tsv is NOT an input here. It was staged and never read --
    // the figure derives its per-platform composition from the per-sample table,
    // because a composition needs the counts and not their summary.
    PLOT_TYPING_QC(TYPING_QC.out.sample_qc
                       .combine(TYPING_QC.out.gene_qc)
                       .combine(RESIDUE_MATRIX.out.sites),
                   script_file('plot_typing_qc.py'))

    // The SECOND panel's summary goes in too. Its tables existed from the day the
    // two-panel fan-out was written and nothing drew them; panel (d) is the whole
    // reason both panels are run, and it was living in a TSV.
    PLOT_ALLELE_FREQ(ch_freq.primary.map   { _m, check, summary, _sm, _pg ->
                         tuple(check, summary) }
                         .combine(
                     ch_freq.secondary.map { _m, _c, summary, _sm, _pg -> summary }),
                     script_file('plot_allele_freq.py'))

    // -- [13] the confound, drawn, under both panels ----------------------------
    PLOT_TYPING_CONFOUND(ch_freq.primary.map   { _m, _c, _s, smp, _pg -> smp }
                             .combine(
                         ch_freq.secondary.map { _m, _c, _s, smp, _pg -> smp }),
                         script_file('plot_typing_confound.py'))

    // -- [14] manifest ---------------------------------------------------------
    def manifest = """{
  "component": "hla.typing",
  "sample_info": "${params.sample_info}",
  "sample_qc_keep": "${params.sample_qc_keep}",
  "cohort_keep": "${params.cohort_keep}",
  "cohort_name": "${params.CohortName}",
  "max_samples": ${params.MaxSamples},
  "reference": "${params.reference}",
  "mhc_window": "${params.MhcChrom}:${params.MhcStart}-${params.MhcEnd}",
  "hlahd": "${params.hlahd_sh}",
  "hlahd_dictionary": "${params.hlahd_dict}",
  "hlahd_gene_split": "${params.hlahd_split}",
  "imgt_alignments": "${params.imgt_alignments}",
  "truth_panel_primary": "${params.hla_truth_jmorp}",
  "truth_panel_secondary": "${params.hla_truth_1kg}",
  "hla_pgroups": "${params.hla_pgroups}",
  "truth_population_1kg": "${params.TruthPopulation}",
  "control_group": "${params.ControlGroup}",
  "reference_loci_primary": ${groovy.json.JsonOutput.toJson(params.RefLociJmorp.split(',') as List)},
  "reference_loci_secondary": ${groovy.json.JsonOutput.toJson(params.RefLoci1kg.split(',') as List)},
  "hlahd_version": "1.7.1",
  "allele_field_depth": ${params.AlleleFieldDepth},
  "min_read_length": ${params.MinReadLength},
  "cutting_rate": ${params.CuttingRate},
  "residue_genes": ${groovy.json.JsonOutput.toJson(params.ResidueGenes.split(',') as List)},
  "min_residue_count": ${params.MinResidueCount}
}"""
    WRITE_RUN_MANIFEST(channel.value(manifest))
}
