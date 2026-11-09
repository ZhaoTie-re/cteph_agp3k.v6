#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The phenotype and covariates for the two models of one cohort and
#           one case set, from one joined table, so the two models can never
#           disagree about who is a case or what a sample's covariates are.
#
#           CASE SET. 'all' keeps every case. 'dvt_yes' / 'dvt_no' keep the cases
#           whose DVT status (case_sets.tsv) is yes / no; the other cases --
#           including every case of unknown status -- leave the analysis. All
#           controls stay in these three sets.
#           'dvt_case_only' drops every control and every case of unknown
#           status, and RECODES PHENO1 within the cases: 1 = DVT yes, 0 = DVT no.
#
#           fixed model  : the cohort's fixed_model sample set (relatives
#                          removed upstream, as in the genome-wide fixed scan).
#                          plink2 reads pheno_covar.fixed.tsv and fixed.keep.
#           random model : every sample of the case set; the GRM in the SAIGE
#                          null model accounts for relatedness. For 'all' the
#                          null model is the genome-wide scan's, reused; for a
#                          DVT set it is fitted from pheno_covar.random.tsv.
#                          SAIGE takes phenotype and covariates from the null
#                          model, and CHECK_NULL proves they match this file.
#
#           PHENO1 is kept in its 0/1 coding (0 = control, 1 = case). plink2
#           reads 0 as MISSING unless it is given --1, so ASSOC_FIXED passes --1
#           and verify.sh checks the case/control counts in every plink2 log
#           against pheno_coding.json.
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import json

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort', required=True)
    p.add_argument('--case-set', required=True, choices=('all', 'dvt_yes', 'dvt_no', 'dvt_case_only'))
    p.add_argument('--case-sets', required=True, help='case_sets.tsv: IID, dvt (yes/no/unknown)')
    p.add_argument('--pheno', required=True)
    p.add_argument('--pheno-col', default='PHENO1')
    p.add_argument('--covar', required=True)
    p.add_argument('--covar-cols', required=True, help='comma-separated, expanded')
    p.add_argument('--fixed-fam', required=True)
    p.add_argument('--out-fixed', required=True)
    p.add_argument('--out-random', required=True)
    p.add_argument('--out-keep', required=True)
    p.add_argument('--out-json', required=True)
    return p.parse_args()


def main():
    a = parse_args()
    cols = [c for c in a.covar_cols.split(',') if c]
    pheno = pd.read_csv(a.pheno, sep='\t', dtype=str)[['IID', a.pheno_col]]
    covar = pd.read_csv(a.covar, sep='\t', dtype=str).rename(columns={'#FID': 'FID'})
    df = pheno.merge(covar[['IID'] + cols], on='IID', how='left', validate='one_to_one')
    if df[cols].isna().any().any():
        raise SystemExit('ERROR: a phenotype sample lacks a covariate value (validate_inputs should have stopped this)')
    df.insert(0, '#FID', df['IID'])                     # the PLINK sets use --double-id

    dvt = pd.read_csv(a.case_sets, sep='\t', dtype=str).set_index('IID')['dvt']
    is_case = df[a.pheno_col] == '1'
    status = df['IID'].map(dvt).fillna('unknown')
    dropped = {}
    n_controls_dropped = 0
    if a.case_set in ('dvt_yes', 'dvt_no'):
        want = a.case_set.split('_')[1]                  # 'yes' / 'no'
        out_cases = is_case & (status != want)
        dropped = status[out_cases].value_counts().to_dict()
        df = df[~out_cases].reset_index(drop=True)
    elif a.case_set == 'dvt_case_only':
        keep = is_case & status.isin(['yes', 'no'])
        dropped = status[is_case & ~keep].value_counts().to_dict()
        n_controls_dropped = int((~is_case).sum())
        df = df[keep].reset_index(drop=True)
        df[a.pheno_col] = (status[keep].values == 'yes').astype(int).astype(str)

    fam = pd.read_csv(a.fixed_fam, sep=r'\s+', header=None, dtype=str)
    fixed = df[df['IID'].isin(set(fam[1]))]

    out_cols = ['#FID', 'IID', a.pheno_col] + cols
    df[out_cols].to_csv(a.out_random, sep='\t', index=False)
    fixed[out_cols].to_csv(a.out_fixed, sep='\t', index=False)
    fixed[['#FID', 'IID']].to_csv(a.out_keep, sep='\t', index=False, header=False)

    def counts(d):
        y = d[a.pheno_col]
        return {'samples': int(len(d)), 'cases': int((y == '1').sum()), 'controls': int((y == '0').sum())}

    coding = ({'0': 'case without DVT', '1': 'case with DVT'} if a.case_set == 'dvt_case_only'
              else {'0': 'control', '1': 'case'})
    info = {'cohort': a.cohort, 'case_set': a.case_set, 'phenotype': a.pheno_col,
            'coding': coding,
            'plink2_flag': '--1',
            'covariates': cols,
            'cases_dropped_by_dvt': {k: int(v) for k, v in dropped.items()},
            'controls_dropped': n_controls_dropped,
            'fixed': counts(fixed), 'random': counts(df)}
    json.dump(info, open(a.out_json, 'w'), indent=1)
    print(f"[prep] {a.cohort} {a.case_set}: dropped cases {info['cases_dropped_by_dvt']} | fixed {info['fixed']} | random {info['random']}")


if __name__ == '__main__':
    main()
