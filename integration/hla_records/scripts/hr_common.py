#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Shared constants and helpers of the table-assembly scripts: cohorts,
#           comparisons and their labels as the decks print them, marker labels,
#           and plain TSV/JSON I/O. The independent re-check (recheck.py) does
#           NOT import this module.
# Component: hla_records
# ---------------------------------------------------------------------------
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

COHORTS = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
COH_LABEL = {'narrow_mainland': 'Narrow mainland', 'intermediate_mainland': 'Intermediate mainland',
             'full_mainland': 'Full mainland'}
COH_CODE = {'narrow_mainland': 'N', 'intermediate_mainland': 'I', 'full_mainland': 'F'}
COMPARISONS = ['all', 'dvt_yes', 'dvt_no', 'dvt_case_only']
# as the association deck names them (CS_LABEL in its build script)
COMP_LABEL = {'all': 'All cases vs. controls', 'dvt_yes': 'With DVT vs. controls',
              'dvt_no': 'Without DVT vs. controls', 'dvt_case_only': 'With vs. without DVT (cases only)'}
COMP_SIDES = {'all': ('All cases', 'Controls'), 'dvt_yes': ('Cases with DVT', 'Controls'),
              'dvt_no': ('Cases without DVT', 'Controls'), 'dvt_case_only': ('Cases with DVT', 'Cases without DVT')}
CLASSES = ['allele', 'residue']
LEAD57, LEAD87, DPB, POS57, POS125 = 'DQB1:57:V', 'DQB1:87:Y', 'DPB1*05:01', 'DQB1:57', 'DQB1:125'
AA = {'A': 'alanine', 'C': 'cysteine', 'D': 'aspartate', 'E': 'glutamate', 'F': 'phenylalanine', 'G': 'glycine',
      'H': 'histidine', 'I': 'isoleucine', 'K': 'lysine', 'L': 'leucine', 'M': 'methionine', 'N': 'asparagine',
      'P': 'proline', 'Q': 'glutamine', 'R': 'arginine', 'S': 'serine', 'T': 'threonine', 'V': 'valine',
      'W': 'tryptophan', 'Y': 'tyrosine'}


def log(*a):
    print(f'[{datetime.now():%H:%M:%S}]', *a, flush=True)


def die(msg):
    sys.exit(f'ERROR: {msg}')


def write_tsv(df, path):
    """Floats at full precision (%.10g), P-like columns as %.6e, '' for missing; never 'nan'."""
    out = df.copy()
    for c in out.columns:
        if pd.api.types.is_float_dtype(out[c]):
            fmt = '%.6e' if (c == 'p' or c.startswith('p_') or c.endswith('_p') or '_p_' in c) else '%.10g'
            out[c] = out[c].map(lambda v: '' if pd.isna(v) else fmt % v)
    out = out.fillna('')
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep='\t', index=False)
    return out


def write_json(obj, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(clean(obj), indent=1) + '\n')


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    return x


def gene_of(marker):
    return marker.split('*')[0] if '*' in marker else marker.split(':')[0]


def marker_label(m):
    """'DQB1:57:V' -> 'DQB1 57 V' (typographic minus), 'DQB1:57' -> 'DQB1 57', alleles unchanged."""
    if '*' in m or ':' not in m:
        return m
    parts = m.split(':')
    parts[1] = parts[1].replace('-', '−')
    return ' '.join(parts)


def gene_rank(genes):
    return {g: i for i, g in enumerate(genes)}
