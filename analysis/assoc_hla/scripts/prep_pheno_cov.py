#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One table carrying the phenotype and every covariate, on exactly the
#           genotyped samples and in the genotype file's order.
#
#           Written here rather than borrowed from assoc_saige so this component
#           has no cross-component runtime dependency: a sibling's refactor must
#           not be able to change what this one fits.
#
#           PHENOTYPE CODING IS NORMALISED, NOT ASSUMED. model_inputs writes 0/1
#           while PLINK's own convention is 1/2, and the two components' READMEs
#           disagree about which. Both are accepted and 0/1 is what leaves here,
#           with the mapping printed.
# Component: assoc_hla
# ---------------------------------------------------------------------------
import argparse
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pheno', required=True, help='model_inputs phenotype/pheno.tsv')
    p.add_argument('--covar', required=True, help='model_inputs covariates/<pc>.sex.tsv')
    p.add_argument('--fam', required=True, help='the genotype .fam; defines the sample set')
    p.add_argument('--pheno-col', default='PHENO1')
    p.add_argument('--id-col', default='IID')
    p.add_argument('--out', default='pheno_covar.tsv')
    return p.parse_args()


def read_tsv(path, id_col):
    d = pd.read_csv(path, sep=r'\s+', dtype=str)
    d.columns = [c.lstrip('#') for c in d.columns]
    if id_col not in d.columns:
        raise SystemExit(f'ABORT: {path} has no {id_col!r} column; found {list(d.columns)[:6]}')
    return d


def main():
    args = parse_args()
    ph = read_tsv(args.pheno, args.id_col)
    cv = read_tsv(args.covar, args.id_col)
    fam = pd.read_csv(args.fam, sep=r'\s+', header=None, dtype=str,
                      names=['FID', 'IID', 'PAT', 'MAT', 'SEX_FAM', 'PHENO_FAM'])

    if args.pheno_col not in ph.columns:
        raise SystemExit(f'ABORT: {args.pheno} has no {args.pheno_col!r}')
    d = fam[['FID', 'IID']].merge(ph[[args.id_col, args.pheno_col]],
                                  left_on='IID', right_on=args.id_col, how='left')
    d = d.drop(columns=[c for c in (args.id_col,) if c != 'IID' and c in d.columns])
    cv2 = cv.drop(columns=[c for c in ('FID',) if c in cv.columns])
    d = d.merge(cv2, left_on='IID', right_on=args.id_col, how='left')

    y = pd.to_numeric(d[args.pheno_col], errors='coerce')
    vals = set(y.dropna().unique())
    if vals <= {1.0, 2.0}:
        d[args.pheno_col] = (y - 1).astype('Int64')
        coding = '1/2 -> 0/1'
    elif vals <= {0.0, 1.0}:
        d[args.pheno_col] = y.astype('Int64')
        coding = '0/1, unchanged'
    else:
        raise SystemExit(f'ABORT: {args.pheno_col} takes {sorted(vals)}; expected 0/1 or 1/2')

    n_missing = int(d[args.pheno_col].isna().sum())
    if n_missing:
        raise SystemExit(f'ABORT: {n_missing} genotyped sample(s) have no phenotype. The '
                         f'cohorts are defined from these files, so this is a join failure, '
                         f'not a real gap -- check for id-type coercion.')
    d.to_csv(args.out, sep='\t', index=False)

    n1 = int((d[args.pheno_col] == 1).sum())
    n0 = int((d[args.pheno_col] == 0).sum())
    print(f'[prep_pheno_cov] {len(d):,} sample(s): {n1:,} case / {n0:,} control '
          f'({coding}) -> {args.out}')
    ncov = [c for c in d.columns if c not in ('FID', 'IID', args.pheno_col)]
    print(f'    {len(ncov)} covariate column(s) available; the model selects from them')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in prep_pheno_cov: {e}', file=sys.stderr)
        sys.exit(1)
