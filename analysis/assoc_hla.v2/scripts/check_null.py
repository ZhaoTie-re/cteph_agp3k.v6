#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Prove that the SAIGE null model a case set's random model will use
#           was fitted on exactly that set: same samples, same case/control label
#           for every sample, covariates SEX + PC1..PC10, binary trait, LOCO.
#
#           For 'all' the null model is the genome-wide scan's, reused; for a
#           DVT set it was fitted by FIT_NULL in this run. The same check runs
#           for both, so every random-model result rests on a checked null.
#           Exits non-zero after writing null_check.json if any check fails.
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import json
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case-set', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--pheno', required=True, help='pheno_covar.random.tsv of the case set')
    p.add_argument('--pheno-col', default='PHENO1')
    p.add_argument('--covar-cols', required=True)
    p.add_argument('--null-samples', required=True)
    p.add_argument('--null-json', required=True)
    p.add_argument('--out', required=True)
    return p.parse_args()


def main():
    a = parse_args()
    ph = pd.read_csv(a.pheno, sep='\t', dtype=str)
    ns = pd.read_csv(a.null_samples, sep='\t', dtype=str)
    nm = json.load(open(a.null_json))
    cols = [c for c in a.covar_cols.split(',') if c]
    checks = []

    def check(name, ok, detail):
        checks.append({'check': name, 'ok': bool(ok), 'detail': detail})

    check('binary trait', nm['trait_type'] == 'binary', nm['trait_type'])
    check('fitted with LOCO', nm['loco'] is True, str(nm['loco']))
    check('sample set equals the case set',
          set(ns.sample_id) == set(ph.IID) and len(ns) == len(ph),
          f'{len(ns)} null samples, {len(ph)} case-set samples, {len(set(ns.sample_id) ^ set(ph.IID))} differ')
    lab = ph.set_index('IID')[a.pheno_col].astype(float)
    ny = ns.set_index('sample_id')['y'].astype(float)
    common = ny.index.intersection(lab.index)
    n_dis = int((ny.loc[common] != lab.loc[common]).sum())
    check('label equals the case set for every sample', n_dis == 0, f'{n_dis} disagree')
    n_case = int((ph[a.pheno_col] == '1').sum())
    check('cases equal the case set', int(nm['cases']) == n_case, f"{nm['cases']} vs {n_case}")
    check('covariates equal the requested list', nm['coefficients'][1:] == cols,
          ','.join(nm['coefficients'][1:]))

    failed = [c for c in checks if not c['ok']]
    json.dump({'case_set': a.case_set, 'cohort': a.cohort, 'n_checks': len(checks),
               'n_failed': len(failed), 'samples': len(ph), 'cases': n_case,
               'variance_ratio': nm['variance_ratio'], 'checks': checks},
              open(a.out, 'w'), indent=1)
    for c in checks:
        print(f"[check_null] {a.case_set} {a.cohort} {'PASS' if c['ok'] else 'FAIL'}  {c['check']}  ({c['detail']})")
    if failed:
        sys.exit(f'ERROR: {len(failed)} null-model checks failed for {a.case_set} {a.cohort}')


if __name__ == '__main__':
    main()
