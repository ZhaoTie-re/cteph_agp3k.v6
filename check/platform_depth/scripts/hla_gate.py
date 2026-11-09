#!/usr/bin/env python3
"""Gate G4 — the HLA re-test machinery on the ORIGINAL matrices and null models reproduces the
published random-model single-marker results (assoc_all.tsv, case set `all`) and the
residue-position omnibus (omnibus.residue.tsv), value for value."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

KEY = ['cohort', 'class', 'marker']
NUM = ['or', 'p', 'p_spa', 'beta', 'se', 'n_case', 'n_control']
OKEY = ['cohort', 'unit']
ONUM = ['chisq', 'p', 'df', 'cases', 'controls']


def diff(ours, pub, key, cols):
    m = ours.merge(pub, on=key, how='outer', suffixes=('', '_pub'), indicator=True)
    out = {'rows_ours': len(ours), 'rows_published': len(pub),
           'only_ours': int((m._merge == 'left_only').sum()), 'only_published': int((m._merge == 'right_only').sum())}
    both = m[m._merge == 'both']
    worst = 0.0
    for c in cols:
        x = pd.to_numeric(both[c], errors='coerce')
        y = pd.to_numeric(both[f'{c}_pub'], errors='coerce')
        nan_mismatch = int((x.isna() != y.isna()).sum())
        ok = x.notna() & y.notna()
        rel = (np.abs(x[ok] - y[ok]) / np.maximum(np.abs(y[ok]), 1e-300)).max() if ok.any() else 0.0
        out[f'max_rel_diff_{c}'] = float(rel)
        out[f'nan_mismatch_{c}'] = nan_mismatch
        worst = max(worst, float(rel), float(nan_mismatch))
    out['worst'] = worst
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--assoc', nargs='+', required=True)
    ap.add_argument('--omnibus', required=True)
    ap.add_argument('--published-assoc', required=True)
    ap.add_argument('--published-omnibus', required=True)
    ap.add_argument('--tol', type=float, default=1e-6)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    ours = pd.concat([pd.read_csv(f, sep='\t', dtype=str) for f in a.assoc], ignore_index=True)
    pub = pd.read_csv(a.published_assoc, sep='\t', dtype=str)
    pub = pub[(pub.case_set == 'all') & (pub.model == 'random')]
    cols = [c for c in NUM if c in ours.columns and c in pub.columns]
    single = diff(ours[KEY + cols], pub[KEY + cols], KEY, cols)

    oo = pd.read_csv(a.omnibus, sep='\t', dtype=str)
    op = pd.read_csv(a.published_omnibus, sep='\t', dtype=str)
    ocols = [c for c in ONUM if c in oo.columns and c in op.columns]
    omni = diff(oo[OKEY + ocols], op[OKEY + ocols], OKEY, ocols)

    ok = (single['only_ours'] == 0 and single['only_published'] == 0 and single['worst'] <= a.tol and
          omni['only_ours'] == 0 and omni['only_published'] == 0 and omni['worst'] <= a.tol)
    rep = {'single_markers': single, 'omnibus': omni, 'compared_columns': {'single': cols, 'omnibus': ocols},
           'pass': bool(ok)}
    Path(a.out).write_text(json.dumps(rep, indent=2) + '\n')
    print(json.dumps(rep, indent=1))
    if not ok:
        sys.exit('GATE G4 FAILED')


if __name__ == '__main__':
    main()
