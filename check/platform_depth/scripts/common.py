"""Shared constants and small helpers for check/platform_depth.

Every script imports this module. Paths are absolute so a step behaves the same on the
login node and inside a SLURM job. Nothing here reads data; it only names it.
"""
from pathlib import Path

PROJECT = Path('/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6')
ROOT = PROJECT / 'check' / 'platform_depth'

COHORTS = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
COHORT_LABEL = {'narrow_mainland': 'Narrow mainland', 'intermediate_mainland': 'Intermediate mainland',
                'full_mainland': 'Full mainland'}

# WGS_Platform labels exactly as the info sheet writes them (note DNBseq vs DNBSeq). The
# slugs are for file names only; tables keep the labels.
PLATFORM_SLUG = {
    'HiSeqX 15x': 'hiseqx',
    'DNBSeq-T7 30x': 't7',
    'NovaSeq 30x': 'novaseq',
    'DNBSeq-G400RS 30x': 'g400rs30',
    'DNBseq-G400RS 15x': 'g400rs15',
}
CONTROL_PLATFORM = 'HiSeqX 15x'
CASE_PLATFORMS = ['DNBSeq-T7 30x', 'NovaSeq 30x', 'DNBSeq-G400RS 30x', 'DNBseq-G400RS 15x']
# The cases down-sampled to the controls' depth (user decision, 2026-10-30).
DEEP_PLATFORMS = ['NovaSeq 30x', 'DNBSeq-G400RS 30x']

# Inputs
INFO_XLSX = PROJECT / 'info' / 'cteph_agp3k.v6.20260507.xlsx'
MODEL_INPUTS = PROJECT / 'wgs.auto.par' / 'results' / '12_model_inputs'
RANDOM_PREFIX = 'cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.random_model'
FIXED_PREFIX = 'cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.maf_ge_threshold'
TUNING_RV = PROJECT / 'tuning.rv' / 'results'
ASSOC_SAIGE = PROJECT / 'analysis' / 'assoc_saige' / 'results'
ASSOC_GENE = PROJECT / 'analysis' / 'assoc_gene' / 'results'
ASSOC_HLA = PROJECT / 'analysis' / 'assoc_hla.v2'
SAIGE_SIGNALS = PROJECT / 'integration' / 'saige_signals' / 'results'
HLA_TYPING = PROJECT / 'hla.typing'
HLA_V2 = PROJECT / 'hla.typing.v2'

FASTA = Path('/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa')
CALLING_INTERVALS = Path('/LARGE1/gr10478/platform/JHRPv6/workspace/pipeline/scripts/interval_list')

# Where per-sample GVCFs of the original calling survive (most were deleted; see docs/METHODS.md).
GVCF_ROOTS = [
    Path('/LARGE0/gr10478/platform/JHRPv5/workspace/pipeline/output/v5'),
    Path('/LARGE1/gr10478/pub/gvcf'),
    Path('/LARGE1/gr10478/workspace/pipeline/output'),
]
# Joint-calling batch lists of the newest (v6) batch; a sample in neither and without a JHRPv5
# output folder is labelled 'v4_or_earlier'. The label only stratifies the re-calling checks.
V6_BATCH_LISTS = [Path('/LARGE1/gr10478/platform/JHRPv6/workspace/pipeline/scripts/'
                       f'{d}/gvcf_lst/gvcf.chr22.lst')
                  for d in ('run_GenomicsDBImport.v6', 'run_GenomicsDBImport.v6.2')]

# Fallback for control CRAM links whose target directory was removed (72 in full mainland).
CRAM_FALLBACK_DIR = Path('/LARGE1/gr10478/pub/WGS/cram/Control_BBJ')
CRAI_FALLBACK_DIR = HLA_TYPING / 'results' / '00.manifest' / 'crai'

# The July audit's genome-wide probe windows (1 Mb each, mid-arm, clear of centromeres).
PROBE_WINDOWS = [(f'chr{c}', 100_000_000, 101_000_000) for c in (1, 2, 3, 4, 5, 6, 7, 8, 10, 11)]

# HaplotypeCaller's read filters, as the cohort ran it (GATK 4.1.x defaults, recorded in the
# cohort VCF header): mapping quality >= 20, base quality >= 10, duplicates excluded.
MIN_MQ = 20
MIN_BQ = 10

# HaplotypeCaller (AssemblyRegionWalker) splits each -L interval into 5,000 bp read shards
# counted from the interval's start; sub-intervals start on that grid so their shards are the
# original run's shards.
SHARD = 5000


def read_tsv(path, **kw):
    import pandas as pd
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, **kw)


def write_tsv(df, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep='\t', index=False)


def split_id(vid):
    """chr<N>:<pos>:<ref>:<alt> -> (chrom, pos, ref, alt)."""
    c, p, r, a = vid.split(':', 3)
    return c, int(p), r, a
