#!/usr/bin/env python3
"""SAIGE step-2 output -> one schema, and gate G1.

Effect and standard error follow the primary, no-SPA P: beta = BETA (= Tstat / var), se = 1/sqrt(var),
so beta / se reproduces p.value.NA. SAIGE's own SE column is not used: on rows where the saddle-point
approximation was applied it equals |BETA| / qnorm(1 - p.value), the one-sided quantile of the
two-sided SPA P.

  collect  stack per-chromosome outputs of one arm x cohort x imputation
  gate     the ORIGINAL arm against the published saige.assoc.tsv (identical printed values)
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

COMPARE = ['BETA', 'SE', 'Tstat', 'var', 'p.value', 'p.value.NA', 'AF_case', 'AF_ctrl', 'N_case', 'N_ctrl',
           'MissingRate', 'AC_Allele2']


def load(files):
    d = pd.concat([pd.read_csv(f, sep=r'\s+', dtype=str) for f in files], ignore_index=True)
    if d.MarkerID.duplicated().any():
        sys.exit('ABORT: a marker appears twice in the SAIGE output')
    return d


def cmd_collect(a):
    d = load(a.saige)
    var = d['var'].astype(float)
    t = d['Tstat'].astype(float)
    out = pd.DataFrame({
        'arm': a.arm, 'cohort': a.cohort, 'subset': a.subset, 'imputation': a.imputation,
        'variant_id': d.MarkerID, 'allele2': d.Allele2,
        'beta': d.BETA.astype(float), 'se': 1 / np.sqrt(var), 'z': t / np.sqrt(var),
        'p': d['p.value.NA'].astype(float), 'p_spa': d['p.value'].astype(float), 'is_spa': d['Is.SPA'],
        'af_case': d.AF_case.astype(float), 'af_ctrl': d.AF_ctrl.astype(float),
        'n_case': d.N_case.astype(int), 'n_ctrl': d.N_ctrl.astype(int),
        'missing_rate': d.MissingRate.astype(float), 'ac_allele2': d.AC_Allele2.astype(float)})
    # beta / se must reproduce the primary P
    pz = 2 * stats.norm.sf(np.abs(out.z))
    bad = np.abs(np.log10(np.maximum(pz, 1e-300)) - np.log10(np.maximum(out.p, 1e-300))) > 1e-3
    if bad.any():
        sys.exit(f'ABORT: {int(bad.sum())} rows where Tstat/sqrt(var) does not give p.value.NA')
    out.to_csv(a.out, sep='\t', index=False)
    print(f'[saige_collect] {a.arm} {a.cohort} {a.subset} {a.imputation}: {len(out)} variants')


def cmd_gate(a):
    ours = load(a.saige).set_index('MarkerID')
    pub = pd.read_csv(a.published, sep='\t', dtype=str, usecols=['MarkerID'] + COMPARE)
    pub = pub[pub.MarkerID.isin(ours.index)].set_index('MarkerID')
    rep = {'cohort': a.cohort, 'variants_ours': len(ours), 'variants_published': len(pub)}
    missing = sorted(set(ours.index) - set(pub.index))
    rep['missing_from_published'] = len(missing)
    diffs = {}
    for c in COMPARE:
        x, y = ours.loc[pub.index, c], pub[c]
        neq = x != y
        if neq.any():
            xf, yf = x[neq].astype(float), y[neq].astype(float)
            rel = (np.abs(xf - yf) / np.maximum(np.abs(yf), 1e-300)).max()
            diffs[c] = {'rows_differing': int(neq.sum()), 'max_rel_diff': float(rel)}
    rep['columns_differing'] = diffs
    rep['max_rel_diff'] = max([v['max_rel_diff'] for v in diffs.values()], default=0.0)
    rep['pass'] = bool(len(missing) == 0 and rep['max_rel_diff'] <= a.tol)
    Path(a.out).write_text(json.dumps(rep, indent=2) + '\n')
    print(json.dumps(rep))
    if not rep['pass']:
        sys.exit('GATE G1 FAILED')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('collect')
    p.add_argument('--saige', nargs='+', required=True)
    p.add_argument('--arm', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--subset', default='pooled')
    p.add_argument('--imputation', required=True)
    p.add_argument('--out', required=True)
    p = sub.add_parser('gate')
    p.add_argument('--saige', nargs='+', required=True)
    p.add_argument('--published', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--tol', type=float, default=1e-6)
    p.add_argument('--out', required=True)
    a = ap.parse_args()
    {'collect': cmd_collect, 'gate': cmd_gate}[a.cmd](a)


if __name__ == '__main__':
    main()
