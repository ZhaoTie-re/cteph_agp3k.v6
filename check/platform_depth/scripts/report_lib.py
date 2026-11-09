"""Shared helpers of the report-stage scripts (ds_summary, platform_report, array_evidence, control_report,
hla_reference_shift, signal_followup, evidence_matrix, figures).

common.py, bedio.py, control_subsets.py and platform_stats.py are staged by cached processes and stay untouched;
what the report scripts share beyond them lives here: the platform, control-source and variant-set names and
labels, the group labels of 13.followup, the signal lists, the replicate arms, the leave-one-out rule, and the
readers of the published per-subset results. Nothing here writes a result.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bedio as B  # noqa: E402
import checks as CK  # noqa: E402
import common as K  # noqa: E402
from control_subsets import CONTROL_SOURCES  # noqa: E402  the info sheet's Project -> source slug

# ── names ────────────────────────────────────────────────────────────────────────────────────────────────────
CASE_SLUGS = [K.PLATFORM_SLUG[p] for p in K.CASE_PLATFORMS]          # t7, novaseq, g400rs30, g400rs15
PLATFORM_ONLY = [f'only_{s}' for s in CASE_SLUGS]                    # one case platform + all controls
PLATFORM_WITHOUT = [f'without_{s}' for s in CASE_SLUGS]              # all case platforms but one + all controls
PLATFORM_SHORT = {'HiSeqX 15x': 'HiSeq X\n(controls)', 'DNBSeq-T7 30x': 'DNBSeq-T7', 'NovaSeq 30x': 'NovaSeq',
                  'DNBSeq-G400RS 30x': 'G400RS\n30x', 'DNBseq-G400RS 15x': 'G400RS\n15x'}
SOURCES = list(CONTROL_SOURCES.values())                             # nagahama, bbj, htlv1, acc, lungca
SOURCE_LABEL = {'nagahama': 'Nagahama', 'bbj': 'BBJ', 'htlv1': 'HTLV-1', 'acc': 'ACC', 'lungca': 'Lung cancer'}
SOURCE_ONLY = [f'ctrlonly_{s}' for s in SOURCES]                     # every case + one source's controls
SOURCE_WITHOUT = [f'ctrlwithout_{s}' for s in SOURCES]               # every case + all sources but one
VARIANT_SETS = ['moderate_high', 'low_moderate_high']                # the genes' variant sets ('stratum' in tables)
SET_SHORT = {'moderate_high': 'M/H', 'low_moderate_high': 'L/M/H'}
COHORT_SHORT = {'narrow_mainland': 'Narrow', 'intermediate_mainland': 'Intermediate', 'full_mainland': 'Full'}
V57 = 'DQB1:57:V'                                                    # the headline HLA residue
POSITION57 = 'DQB1:57'                                               # its residue position (omnibus unit)
FULL = 'full_mainland'                                                # the cohort that holds the other two

# ── group labels of 13.followup (written by signal_followup.py, read by evidence_matrix.py and figures.py):
#    'case | <WGS_Platform>[ | <new|old> study]', 'control | <source>', the pooled groups, the reference panels
CASE, CONTROL, REFERENCE = 'case | ', 'control | ', 'reference | '
CASE_ALL, CONTROL_ALL = CASE + 'all', CONTROL + 'all'
REF_TOMMO, REF_JPT = REFERENCE + 'ToMMo 61KJPN', REFERENCE + '1KG JPT'
STUDY = {'新研究': 'new', '旧研究': 'old'}                            # info sheet Registration_Study


# ── signals (02.targets/signal_variants.tsv) ────────────────────────────────────────────────────────────────
def read_signals(results):
    return pd.read_csv(Path(results) / '02.targets' / 'signal_variants.tsv', sep='\t', dtype=str)


def gwas_leads(sig, alt=False):
    """the locus leads (and the alternative leads), one row per variant, in table order"""
    roles = ['Lead', 'Alt lead'] if alt else ['Lead']
    return sig[(sig.kind == 'gwas') & sig.role.isin(roles)].drop_duplicates('variant_id')


def saige_signals(sig):
    """(unit, marker) of every signal tested by SAIGE: the GWAS leads, then DQB1:57:V"""
    return gwas_leads(sig)[['unit', 'variant_id']].values.tolist() + [[V57, V57]]


def gene_sets(sig):
    """(gene, variant set) of the Tier 1-2 genes, sorted"""
    return sig[sig.kind == 'gene'][['unit', 'stratum']].drop_duplicates().sort_values(['unit', 'stratum']).values.tolist()


# ── replicate arms ───────────────────────────────────────────────────────────────────────────────────────────
def replicate_arms(columns, prefix):
    """r01, r02, ... from the columns <prefix>r01, <prefix>r02, ... of a ds.* table"""
    return sorted(c[len(prefix):] for c in columns
                  if c.startswith(prefix) and c[len(prefix):][:1] == 'r' and c[len(prefix) + 1:].isdigit())


# ── z and the leave-one-out rule (METHODS §11; thresholds in checks.py) ────────────────────────────────────
def z_from_p(p):
    """|z| of a two-sided P"""
    return norm.isf(np.asarray(p, float) / 2)


def loo_flags(z_obs, z_exp):
    """elementwise: testable (|z expected| >= 2.8) and |z observed| < 1.96 or the opposite sign"""
    zo, ze = np.asarray(z_obs, float), np.asarray(z_exp, float)
    return (np.abs(ze) >= CK.LOO_Z_EXPECTED) & ((np.abs(zo) < CK.LOO_Z_OBSERVED) | (np.sign(zo) != np.sign(ze)))


def loo_reading(rows):
    """rows: (label, z_obs, z_exp); returns the reading and the flagged (or else the weakest) left-out groups"""
    t = [(lab, zo, ze) for lab, zo, ze in rows if np.isfinite(ze) and abs(ze) >= CK.LOO_Z_EXPECTED and np.isfinite(zo)]
    if not t:
        return 'not testable', ''
    bad = [(lab, zo, ze) for lab, zo, ze in t if abs(zo) < CK.LOO_Z_OBSERVED or np.sign(zo) != np.sign(ze)]
    w = bad or sorted(t, key=lambda r: abs(r[1]) - abs(r[2]))[:1]
    txt = '; '.join(f'-{lab}: z {zo:.2f} (expected {ze:.2f})' for lab, zo, ze in w)
    return ('flag' if bad else 'consistent'), txt


# ── published per-subset results (08.subsets, 09.platform_assoc, 10.platform) ──────────────────────────────
def platform_table(results, cohort, name):
    return Path(results) / '10.platform' / cohort / f'platform.{name}.tsv'


def cohorts_with(results, cohort_file):
    """the cohorts whose file cohort_file(results, cohort) exists, in cohort order"""
    return [c for c in K.COHORTS if cohort_file(results, c).exists()]


def subset_counts_file(results, cohort, subset):
    """08.subsets/<cohort>/{subsets|control_sources/subsets}/subset_counts.tsv, the one that holds <subset>"""
    d = Path(results) / '08.subsets' / cohort
    return d / ('control_sources/subsets' if subset.startswith('ctrl') else 'subsets') / 'subset_counts.tsv'


def read_saige_subsets(results, cohort, subsets):
    """{subset: 09.platform_assoc/<cohort>/<subset>.saige.tsv, one row per variant_id} for the subsets present"""
    d = {}
    for sub in subsets:
        f = Path(results) / '09.platform_assoc' / cohort / f'{sub}.saige.tsv'
        if f.exists():
            d[sub] = pd.read_csv(f, sep='\t', dtype={'variant_id': str}).drop_duplicates('variant_id').set_index('variant_id')
    return d


def gene_subset_file(results, cohort, subset, vset):
    return Path(results) / '09.platform_assoc' / cohort / f'gene.{subset}.{vset}.tsv'


def omnibus_file(results, subset):
    return Path(results) / '09.platform_assoc' / '_omnibus' / f'omnibus.{subset}.tsv'


def omnibus_wide(files):
    """[(subset, omnibus.<subset>.tsv)] -> (the long table, P per cohort x unit with one p_<subset> column each)"""
    om = []
    for sub, f in files:
        x = pd.read_csv(f, sep='\t', dtype=str)
        x['subset'] = sub
        om.append(x)
    if not om:
        return None, None
    om = pd.concat(om)
    om['p'] = om.p.astype(float)
    w = om.pivot_table(index=['cohort', 'unit'], columns='subset', values='p', aggfunc='first')
    w.columns = [f'p_{s}' for s in w.columns]
    return om, w.reset_index()


# ── genotypes ────────────────────────────────────────────────────────────────────────────────────────────────
def alt_dosage(prefix):
    """IIDs, variant IDs and ALT copies (float, NaN = missing; variants x samples) of a PLINK 1 fileset, read by
    bedio.py; the ALT is the last field of the CHROM:POS:REF:ALT variant ID"""
    fam, bim = B.read_fam(prefix), B.read_bim(prefix)
    codes = B.unpack(B.read_rows(prefix, list(range(len(bim))), len(fam)), len(fam))
    col5 = np.select([codes == B.HOM5, codes == B.HET, codes == B.HOM6], [2.0, 1.0, 0.0], np.nan)
    col5_is_alt = np.array([f[4] == f[1].split(':')[3] for f in bim])
    dose = np.where(col5_is_alt[:, None], col5, 2 - col5)
    return np.array(fam, dtype=object), np.array([f[1] for f in bim], dtype=object), dose
