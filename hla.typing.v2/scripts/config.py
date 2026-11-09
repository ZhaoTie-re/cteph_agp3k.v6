"""Single source of paths, versions and thresholds for hla.typing.v2.

Every step imports this module. Paths are absolute so a step behaves the same on
the login node and inside a SLURM job. Thresholds are written here BEFORE the step
that uses them is run; their fixing date is recorded beside each group.
"""
from pathlib import Path

# --------------------------------------------------------------------------
# Locations
# --------------------------------------------------------------------------
PROJECT = Path('/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6')
ROOT = PROJECT / 'hla.typing.v2'
RESULTS = ROOT / 'results'
D_INPUTS = RESULTS / '00.inputs'
D_CALLS = RESULTS / '01.calls'
D_EVIDENCE = RESULTS / '02.evidence'
D_PANELS = RESULTS / '03.panels'
D_CONF = RESULTS / '04.confidence'
D_EXPORTS = RESULTS / '05.exports'
D_ASSOC = RESULTS / '06.assoc_inputs'
D_FIG = RESULTS / 'figures'

# The earlier typing run (HLA-HD over the same 3,101 samples, 2026-08-25). Its
# per-sample output directories are symlinks into that component's work/ tree.
OLD = PROJECT / 'hla.typing'
OLD_TYPING = OLD / 'results' / '02.typing'
OLD_READS = OLD / 'results' / '01.reads'
OLD_MANIFEST = OLD / 'results' / '00.manifest' / 'sample_manifest.tsv'
OLD_CONTIGS = OLD / 'results' / '00.manifest' / 'contig_list.txt'

# The cohort definition (ancestry filter; 3,101 = 2,662 AGP3K + 439 PH).
COHORT_KEEP = PROJECT / 'PopGMM_output' / 'full_mainland.fid_iid.txt'
CONTROL_GROUP = 'AGP3K'

# --------------------------------------------------------------------------
# Software and references used by the typing run (recorded, not re-run)
# --------------------------------------------------------------------------
HLAHD_DIR = Path('/LARGE0/gr10478/bin/hlahd.1.7.1')
HLAHD_FREQ = HLAHD_DIR / 'freq_data'
HLAHD_SPLIT = Path('/LARGE0/gr10478/b37974/software/hlahd.dictionary/HLA_gene.split.3.50.0.txt')
HLAHD_DICT = Path('/LARGE0/gr10478/b37974/software/hlahd.dictionary/IMGT-HLA-3.64.0_N150')
IMGT_ALIGN = Path('/LARGE0/gr10478/b37974/software/hlahd.dictionary/IMGT-alignments-3.64.0')
P_GROUPS = PROJECT.parent / 'jMorp_HLA_types' / 'hla_nom_p.txt'          # IPD-IMGT 3.64.0
# Every allele name in every IPD-IMGT/HLA release up to 3.64.0 (GitHub ANHIG/IMGTHLA,
# tag 3640), pinned: dates when each allele was first named, so that an allele the
# panel's older dictionary could not contain is not counted against our calls.
IMGT_HISTORY = ROOT / 'assets' / 'imgt' / 'Allelelist_history.txt'

# --------------------------------------------------------------------------
# External panels
# --------------------------------------------------------------------------
TOMMO = PROJECT.parent / 'jMorp_HLA_types' / 'HLA_allele_frequencies_61K.txt'
TOMMO_MD5 = PROJECT.parent / 'jMorp_HLA_types' / 'MD5SUM'
KG = PROJECT.parent / '1KG_HLA_types' / '20181129_HLA_types_full_1000_Genomes_Project_panel.txt'
KG_MD5 = PROJECT.parent / '1KG_HLA_types' / 'MD5SUM'
KG_POPULATION = 'JPT'
TOMMO_LOCI = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1', 'DRB3', 'DRB4', 'E', 'F', 'G']
KG_LOCI = ['A', 'B', 'C', 'DQB1', 'DRB1']
# ToMMo loci whose allele frequencies are NOT compared (the row stays in locus_concordance.tsv with
# frequency_comparable = 0 and this reason). Their absent-allele check still runs.
TOMMO_FREQ_NOT_COMPARABLE = {
    'DRB3': 'not comparable for frequencies: the panel assigns a DRB3 allele to every chromosome, including '
            'the haplotypes that carry no DRB3',
    'DRB4': 'not compared for frequencies: presence follows DRB1 (DR53); the panel\'s non-null DRB4 share '
            '(45.7 %) does not match the DR53 share of its own DRB1 (35.3 %); and among carriers most alleles '
            'are DRB4*01:01 (92 % here, 99 % in the panel), so r carries little information',
}

# --------------------------------------------------------------------------
# Interpreter (absolute: `source activate` does not put the env first on PATH)
# --------------------------------------------------------------------------
PYTHON = '/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3'
PLINK2 = '/home/b/b37974/plink2'    # only for the deck's genotype comparison (slide 19)

# --------------------------------------------------------------------------
# Genes
# --------------------------------------------------------------------------
# HLA-HD types the 33 genes of its gene-split file. 19 are functional genes and
# are reported; the other 14 are pseudogenes, kept in the raw tables only.
REPORTED_GENES = ['A', 'B', 'C', 'E', 'F', 'G',
                  'DRA', 'DRB1', 'DRB3', 'DRB4', 'DRB5',
                  'DQA1', 'DQB1', 'DPA1', 'DPB1',
                  'DMA', 'DMB', 'DOA', 'DOB']
# Present on every haplotype: a missing call is a failure, a single allele is a
# homozygote.
OBLIGATE_GENES = [g for g in REPORTED_GENES if g not in ('DRB3', 'DRB4', 'DRB5')]
# Present only on some DR haplotypes; copy number follows DRB1 (DR51/52/53).
PARALOGUES = ['DRB3', 'DRB4', 'DRB5']
CLASS_I = ['A', 'B', 'C', 'E', 'F', 'G']
# States that carry no genotype.
NO_CALL_STATES = ['failed', 'unresolved', 'missing']

# DRB1 allele group (first field) -> the paralogue its haplotype carries.
# DR51: DRB1*15, *16 -> DRB5. DR52: DRB1*03, *11, *12, *13, *14 -> DRB3.
# DR53: DRB1*04, *07, *09 -> DRB4. DR1, DR8, DR10 (DRB1*01, *08, *10): none.
DR_GROUP_PARALOGUE = {
    '15': 'DRB5', '16': 'DRB5',
    '03': 'DRB3', '11': 'DRB3', '12': 'DRB3', '13': 'DRB3', '14': 'DRB3',
    '04': 'DRB4', '07': 'DRB4', '09': 'DRB4',
    '01': None, '08': None, '10': None,
}

# --------------------------------------------------------------------------
# Thresholds — per call (fixed 2026-10-19, before 02.evidence was run)
# --------------------------------------------------------------------------
MIN_ALLELE_READS = 10        # an allele supported by fewer reads is `low_support`
MIN_MINOR_SHARE = 0.15       # minor allele's share of a heterozygote's UNIQUE reads; below = `imbalanced`
MIN_UNIQUE_READS = 10        # balance is assessed only when the two alleles have this many unique reads
# A tie between allele pairs broken by HLA-HD's frequency prior counts as resolved
# when the chosen pair's prior weight is at least this many times the best other
# pair's (weight = product of the two alleles' counts in HLA-HD's freq_data).
PRIOR_DOMINANCE = 100

# --------------------------------------------------------------------------
# Thresholds — per sample (fixed 2026-10-19, before 02.evidence was run)
# --------------------------------------------------------------------------
OUTLIER_READ_PAIRS_PCTL = 1.0      # read pairs below this cohort percentile
OUTLIER_MIN_UNRESOLVED = 3         # this many unresolved reported genes or more
# (any failed obligate gene also makes a sample an outlier)

# --------------------------------------------------------------------------
# Thresholds — per locus
# --------------------------------------------------------------------------
# Fixed 2026-10-19 before 03.panels was run, and kept:
TIER_MAX_UNCONFIRMED = 0.02        # share of control chromosomes on an allele absent from ToMMo
                                   # (counted over alleles ToMMo's dictionary could contain;
                                   # see nomenclature.first_release and panels.tommo_vintage)
TIER_MIN_HIGH_CONF = 0.95          # share of calls with high confidence
# Revised 2026-10-19 after the first run of 03.panels (docs/METHODS.md §6, revision record):
#  - the absent-allele share counts only alleles ToMMo's dictionary could contain (first
#    named by IPD-IMGT no later than ToMMo's newest allele, release 3.21.0); the 880 of 1,355
#    absent chromosomes on later alleles are reported separately. The 2 % threshold is unchanged.
#  - frequency discordance against ToMMo is REPORTED, NOT TIERED. Among the alleles both
#    sources carry it did not shrink on high-confidence calls, and 1KG JPT sided with our
#    controls on 24 of the 33 alleles shifted by >= 1 point: it mixes cohort differences with typing.
#    The value that was to be its threshold is kept for the figures' reference line.
#  - homozygosity is judged as the share of heterozygotes reported as homozygotes
#    (dropout = (observed - expected) / (1 - expected)); a locus fails when the lower
#    95 % bound exceeds TIER_MAX_DROPOUT, the 3 % tolerance the frequency criterion had.
#    The first rule ("the CI of observed/expected includes 1") tolerated no excess at all.
TIER_MAX_EXCESS_D = 0.03           # reference line only
TIER_MAX_DROPOUT = 0.03
N_SIM = 1000                       # sampling draws for the expected discordance
SEED = 20261019

# --------------------------------------------------------------------------
# Association inputs (step 06 writes them; no test is run here; floors chosen 2026-10-20)
# --------------------------------------------------------------------------
# The association runs in three nested cohorts, each with its own phenotype file; a marker's call
# rate and frequency are measured within the cohort it is tested in, so each cohort gets its own matrices.
COHORTS = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
MODEL_INPUTS = PROJECT / 'wgs.auto.par' / 'results' / '12_model_inputs'   # <cohort>/phenotype/pheno.tsv
ASSOC_GENES = OBLIGATE_GENES   # DRB3, DRB4 and DRB5 left out: their presence follows DRB1
MIN_CALL_RATE = 0.95           # share with a genotype at the marker, required in cases AND in controls (the
                               # sequencing platform and depth differ between them); applied first
MIN_MAF = 0.005                # copies of the less common form over the copies typed; applied second
