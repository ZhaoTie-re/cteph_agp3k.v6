#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The one vocabulary every assoc_hla.v2 script and figure uses.
#
#           Names, orders and colours that mean something in this component are
#           defined here once, so a table column, a figure legend and a README
#           sentence never call the same thing two names. The words follow the
#           hla.typing.v2 deck: a marker is one column of a 0/1/2 matrix; an
#           "allele" column counts copies of one 2-field allele, a "residue"
#           column counts copies carrying one amino acid at one IMGT position;
#           cohorts are "narrow / intermediate / full mainland".
#
#           COLOUR BY MEANING, from the shared palette: the two models are
#           SERIES[0] (fixed) and SERIES[1] (random); cohorts, an ordered
#           nesting, take COHORT_RAMP; thresholds are ACCENT (GEC, primary),
#           REFERENCE (max-T) and INK_SOFT (conventional Bonferroni).
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
import plot_style as S                 # noqa: E402

# ── design ───────────────────────────────────────────────────────────────────
COHORTS = ('narrow_mainland', 'intermediate_mainland', 'full_mainland')
COHORT_LABEL = {'narrow_mainland': 'narrow mainland',
                'intermediate_mainland': 'intermediate mainland',
                'full_mainland': 'full mainland'}
COHORT_SHORT = {'narrow_mainland': 'narrow', 'intermediate_mainland': 'intermediate',
                'full_mainland': 'full'}
COHORT_COLOR = dict(zip(COHORTS, S.COHORT_RAMP))

# Case definitions. 'all', 'dvt_yes' and 'dvt_no' compare their cases with ALL
# the cohort's controls; a case whose DVT status is unknown enters 'all' only.
# 'dvt_case_only' has no controls: it compares cases with DVT (coded 1) against
# cases without DVT (coded 0), the direct test of a difference between the two.
# In its tables the "case" columns are DVT yes and the "control" columns DVT no.
CASE_SETS = ('all', 'dvt_yes', 'dvt_no', 'dvt_case_only')
CASE_SET_LABEL = {'all': 'CTEPH (all cases) vs controls',
                  'dvt_yes': 'CTEPH with DVT vs controls',
                  'dvt_no': 'CTEPH without DVT vs controls',
                  'dvt_case_only': 'CTEPH with DVT vs CTEPH without DVT (cases only)'}
CASE_SET_SHORT = {'all': 'all cases', 'dvt_yes': 'DVT yes', 'dvt_no': 'DVT no',
                  'dvt_case_only': 'DVT yes vs no, cases only'}
CASE_SET_COLOR = {'all': S.REFERENCE, 'dvt_yes': S.SERIES[1], 'dvt_no': S.SERIES[0],
                  'dvt_case_only': S.INK}

CLASSES = ('allele', 'residue')
CLASS_LABEL = {'allele': 'allele', 'residue': 'amino-acid residue'}
CLASS_PLURAL = {'allele': 'alleles', 'residue': 'amino-acid residues'}

MODELS = ('fixed', 'random')
MODEL_LABEL = {'fixed': 'fixed model (logistic regression)',
               'random': 'random model (logistic mixed model)'}
MODEL_SHORT = {'fixed': 'fixed', 'random': 'random'}
MODEL_ENGINE = {'fixed': 'plink2 --glm', 'random': 'SAIGE step 2'}
MODEL_COLOR = {'fixed': S.SERIES[0], 'random': S.SERIES[1]}
MODEL_MARKER = {'fixed': 'o', 'random': 's'}

# ── significance: three family-wise thresholds per (case set, cohort, class) ─
# All at alpha = 0.05 over the class's tests (a two-residue mirror pair is one
# test), fixed from the genotypes before any P exists (effective_tests.py).
THRESHOLDS = ('gec', 'maxt', 'bonf')
THRESHOLD_LABEL = {'gec': 'LD-adjusted (GEC)', 'maxt': 'Exact (max-T)',
                   'bonf': 'Conventional Bonferroni'}
THRESHOLD_ROLE = {'gec': 'primary', 'maxt': 'validation', 'bonf': 'sensitivity'}
THRESHOLD_COLOR = {'gec': S.ACCENT, 'maxt': S.REFERENCE, 'bonf': S.INK_SOFT}
THRESHOLD_LS = {'gec': '-', 'maxt': '--', 'bonf': ':'}

# Cross-cohort consistency: in how many of the three nested cohorts a marker
# reaches a threshold in the random (primary) model, with the same direction of
# effect. A sample-selection sensitivity analysis, not replication. 1/3 is
# called "Single cohort", as in the gene-based analysis.
PRIMARY_MODEL = 'random'


# Its colour: one hue, light -> dark with n, and 0/3 a recessive grey
# (dataviz validator: adjacent steps >= 16 dE under CVD, >= 19 normal vision).
# Mark size grows with n as a second encoding, and signal leads are labelled.
CONS_COLOR = {0: '#CDD2D7', 1: '#E3896B', 2: '#B53F3B', 3: '#5E1A18'}
CONS_LABEL = {3: '3/3 cohorts', 2: '2/3 cohorts', 1: '1/3 (Single cohort)', 0: 'not significant'}
CONS_WIDTH = {0: 0.6, 1: 1.3, 2: 1.6, 3: 1.9}        # range-bar line width (points)
CONS_SIZE = {0: 4, 1: 12, 2: 16, 3: 20}              # scatter marker area


def consistency(n, of=3):
    """'3/3', '2/3', '1/3', '0/3'."""
    return f'{int(n)}/{of}'


def consistency_words(n, of=3):
    """Text form: 'significant in 3/3 cohorts', 'Single cohort', 'not significant'."""
    n = int(n)
    if n == 0:
        return 'not significant'
    if n == 1:
        return 'Single cohort'
    return f'significant in {n}/{of} cohorts'


# ── genes: chromosomal order on GRCh38 (set at run time from the coordinates) ─
GENES = ('A', 'B', 'C', 'E', 'F', 'G', 'DRA', 'DRB1', 'DQA1', 'DQB1',
         'DPA1', 'DPB1', 'DMA', 'DMB', 'DOA', 'DOB')


def gene_symbol(gene):
    """IPD-IMGT/HLA locus name -> HGNC symbol ('DRB1' -> 'HLA-DRB1')."""
    return f'HLA-{gene}'


# ── marker names ─────────────────────────────────────────────────────────────
# allele : 'A*02:01'        gene 'A', no position
# residue: 'DRB1:-27:L'     gene 'DRB1', IMGT position -27, residue 'L'
_RESIDUE = re.compile(r'^([A-Z0-9]+):(-?\d+):([A-Z*.-])$')
_ALLELE = re.compile(r'^([A-Z0-9]+)\*(\d+):(\d+)$')


def parse_marker(name, cls):
    """(gene, position or None, residue or None) for one matrix column name."""
    if cls == 'residue':
        m = _RESIDUE.match(name)
        if not m:
            raise ValueError(f'not a residue column name: {name!r}')
        return m.group(1), int(m.group(2)), m.group(3)
    m = _ALLELE.match(name)
    if not m:
        raise ValueError(f'not a 2-field allele column name: {name!r}')
    return m.group(1), None, None


def sort_key(name, cls):
    """Order of markers within one gene: residues by IMGT position then residue
    letter; alleles by the two numeric fields."""
    if cls == 'residue':
        _, pos, aa = parse_marker(name, cls)
        return (pos, aa)
    m = _ALLELE.match(name)
    return (int(m.group(2)), int(m.group(3)))


def marker_label(name, cls):
    """Display form: alleles keep their name; residues read 'DRB1 position 13 (S)'."""
    if cls == 'residue':
        gene, pos, aa = parse_marker(name, cls)
        return f'{gene} {pos} {aa}'
    return name
