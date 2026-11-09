#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Check, per cohort, that every input this component reads is the one
#           the design assumes -- before any test is run.
#
#           The matrices are hla.typing.v2's step-06 export, verified there end
#           to end; this step proves the copy read here is that export (md5
#           against summary.json) and that it lines up with the phenotype file
#           row for row. The SAIGE null model is reused from the genome-wide
#           scan, so its samples, case labels and covariates are compared with
#           this component's own phenotype and covariate files: SAIGE step 2
#           matches samples by ID and silently drops what it cannot match.
#
#           The DVT status every case set is cut from must cover every case.
#
#           Every check is written to validation.json with its outcome; the
#           script exits non-zero if any check failed, AFTER writing the file,
#           so the failure is readable on disk.
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import hashlib
import json
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort', required=True)
    p.add_argument('--summary', required=True, help='06.assoc_inputs/summary.json')
    p.add_argument('--markers', required=True, help='06.assoc_inputs/markers.tsv')
    p.add_argument('--allele', required=True)
    p.add_argument('--residue', required=True)
    p.add_argument('--pheno', required=True)
    p.add_argument('--pheno-col', default='PHENO1')
    p.add_argument('--covar', required=True)
    p.add_argument('--covar-cols', required=True, help='comma-separated, expanded')
    p.add_argument('--fixed-fam', required=True)
    p.add_argument('--null-samples', required=True)
    p.add_argument('--null-json', required=True)
    p.add_argument('--case-sets', required=True, help='case_sets.tsv: IID, dvt (yes/no/unknown)')
    p.add_argument('--out', required=True)
    return p.parse_args()


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    a = parse_args()
    checks = []

    def check(name, ok, detail):
        checks.append({'check': name, 'ok': bool(ok), 'detail': detail})

    summ = json.load(open(a.summary))['cohorts'][a.cohort]
    markers = pd.read_csv(a.markers, sep='\t', dtype=str, keep_default_na=False)
    pheno = pd.read_csv(a.pheno, sep='\t', dtype=str)
    covar = pd.read_csv(a.covar, sep='\t', dtype=str)
    covar = covar.rename(columns={'#FID': 'FID'})
    covar_cols = [c for c in a.covar_cols.split(',') if c]

    # ── 1. the matrices are the export the typing component verified ─────────
    want_n = {'allele': int(summ['alleles_tested']), 'residue': int(summ['residues_tested'])}
    for cls, path in (('allele', a.allele), ('residue', a.residue)):
        got = md5(path)
        want = summ['md5'][f'{cls}_dosage.tsv']
        check(f'{cls}: md5 equals summary.json', got == want, f'{got} vs {want}')

        m = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
        cols = list(m.columns[1:])
        tested = markers.loc[(markers.cohort == a.cohort) & (markers['class'] == cls)
                             & (markers.status == 'tested'), 'marker'].tolist()
        check(f'{cls}: column count equals summary.json', len(cols) == want_n[cls],
              f'{len(cols)} vs {want_n[cls]}')
        check(f'{cls}: columns equal the markers.tsv "tested" rows', sorted(cols) == sorted(tested),
              f'{len(set(cols) ^ set(tested))} names differ')
        check(f'{cls}: column names unique', len(set(cols)) == len(cols), f'{len(cols)} columns')
        check(f'{cls}: samples equal the phenotype file, in order',
              m['sample_id'].tolist() == pheno['IID'].tolist(),
              f'{len(m)} matrix rows, {len(pheno)} phenotype rows')
        vals = set(pd.unique(m.iloc[:, 1:].values.ravel()))
        check(f'{cls}: values in {{blank, 0, 1, 2}}', vals <= {'', '0', '1', '2'},
              'values: ' + ', '.join(sorted(repr(v) for v in vals)))

    # ── 2. phenotype ─────────────────────────────────────────────────────────
    y = pheno[a.pheno_col]
    check('phenotype: coded 0/1, no missing', set(y) <= {'0', '1'},
          'values: ' + ', '.join(sorted(set(y))))
    n_case, n_ctrl = int((y == '1').sum()), int((y == '0').sum())
    check('phenotype: samples equal summary.json', len(pheno) == int(summ['samples']),
          f'{len(pheno)} vs {summ["samples"]}')
    check('phenotype: cases equal summary.json', n_case == int(summ['cases']),
          f'{n_case} vs {summ["cases"]}')
    check('phenotype: IID unique', pheno['IID'].is_unique, f'{len(pheno)} rows')

    # ── 3. covariates: all of them, for every sample ─────────────────────────
    miss_cols = [c for c in covar_cols if c not in covar.columns]
    check('covariates: every requested column present', not miss_cols,
          'missing: ' + ', '.join(miss_cols) if miss_cols else ', '.join(covar_cols))
    cv = covar.set_index('IID').reindex(pheno['IID'])
    n_nosample = int(cv.isna().all(axis=1).sum())
    check('covariates: every phenotype sample has a row', n_nosample == 0, f'{n_nosample} without a row')
    if not miss_cols:
        num = cv[covar_cols].apply(pd.to_numeric, errors='coerce')
        n_bad = int(num.isna().any(axis=1).sum())
        check('covariates: all values numeric and present', n_bad == 0, f'{n_bad} samples incomplete')
        sexes = sorted(num['SEX'].dropna().unique().tolist()) if 'SEX' in num else []
        check('covariates: SEX takes two values', len(sexes) == 2, f'values {sexes}')

    # ── 4. fixed-model sample set: relatives removed, a subset of the cohort ─
    fam = pd.read_csv(a.fixed_fam, sep=r'\s+', header=None, dtype=str)
    fam_iid = set(fam[1])
    outside = len(fam_iid - set(pheno['IID']))
    check('fixed set: every .fam sample is in the phenotype file', outside == 0, f'{outside} outside')
    check('fixed set: smaller than the cohort (relatives removed)', len(fam_iid) < len(pheno),
          f'{len(fam_iid)} of {len(pheno)}')

    # ── 5. the reused null model was fitted on exactly these samples ─────────
    nm = json.load(open(a.null_json))
    ns = pd.read_csv(a.null_samples, sep='\t', dtype=str)
    check('null model: binary trait', nm['trait_type'] == 'binary', nm['trait_type'])
    check('null model: fitted with LOCO', nm['loco'] is True, str(nm['loco']))
    check('null model: sample set equals the phenotype file',
          set(ns['sample_id']) == set(pheno['IID']) and len(ns) == len(pheno),
          f'{len(ns)} null samples, {len(set(ns.sample_id) ^ set(pheno.IID))} differ')
    lab = pheno.set_index('IID')[a.pheno_col].astype(float)
    ny = ns.set_index('sample_id')['y'].astype(float)
    common = ny.index.intersection(lab.index)
    n_disagree = int((ny.loc[common] != lab.loc[common]).sum())
    check('null model: case/control label equals the phenotype file for every sample',
          n_disagree == 0, f'{n_disagree} disagree')
    check('null model: cases equal the phenotype file', int(nm['cases']) == n_case,
          f'{nm["cases"]} vs {n_case}')
    check('null model: covariates equal the requested list',
          nm['coefficients'][1:] == covar_cols,
          'null: ' + ','.join(nm['coefficients'][1:]))

    # ── 6. DVT status: every case is in the sample sheet ─────────────────────
    cs = pd.read_csv(a.case_sets, sep='\t', dtype=str).set_index('IID')['dvt']
    case_ids = pheno.loc[y == '1', 'IID']
    n_absent = int((~case_ids.isin(cs.index)).sum())
    check('DVT: every case has a row in the sample sheet', n_absent == 0, f'{n_absent} cases absent')
    dvt_counts = case_ids.map(cs).fillna('unknown').value_counts().to_dict()
    check('DVT: controls carry no DVT status',
          not pheno.loc[y == '0', 'IID'].map(cs).isin(['yes', 'no']).any(), 'controls with yes/no')

    failed = [c for c in checks if not c['ok']]
    out = {'cohort': a.cohort, 'n_checks': len(checks), 'n_failed': len(failed),
           'samples': len(pheno), 'cases': n_case, 'controls': n_ctrl,
           'fixed_samples': len(fam_iid), 'variance_ratio': nm['variance_ratio'],
           'case_dvt': {k: int(v) for k, v in dvt_counts.items()},
           'checks': checks}
    json.dump(out, open(a.out, 'w'), indent=1)
    for c in checks:
        print(f"[validate] {a.cohort} {'PASS' if c['ok'] else 'FAIL'}  {c['check']}  ({c['detail']})")
    if failed:
        sys.exit(f'ERROR: {len(failed)} of {len(checks)} input checks failed for {a.cohort}')


if __name__ == '__main__':
    main()
