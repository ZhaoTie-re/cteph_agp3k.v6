#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One table for SAIGE step 1: phenotype and every covariate, on the
#           samples that are actually in the genotypes.
#
#           SAIGE takes a SINGLE --phenoFile carrying both the phenotype and the
#           covariate columns, while the upstream keeps them apart (one phenotype
#           file, one covariate file per PC space). This joins them on the sample
#           ID and writes every covariate column through, because the covariate
#           SUBSET is chosen per model by --covarColList — the same table serves
#           a model fitted with PCs and one fitted without.
#
#           --keep-fam restricts to the genotype samples. SAIGE would intersect
#           anyway, but doing it here means the row count in this file is the N
#           the null model will actually fit, so a silent sample loss shows up
#           now rather than as a puzzling N in the log hours later.
# Component: assoc_saige
# Used by : assoc_saige.nf  process PREP_PHENO_COV
# ---------------------------------------------------------------------------
import argparse
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description='Phenotype + covariates -> one SAIGE step-1 table.')
    p.add_argument('--pheno', required=True)
    p.add_argument('--cov', required=True)
    p.add_argument('--pheno-col', required=True)
    p.add_argument('--id-col', default='IID')
    p.add_argument('--keep-fam', help='plink .fam; restrict to these samples.')
    p.add_argument('--out', required=True)
    return p.parse_args()


def read_table(path, what):
    d = pd.read_csv(path, sep=r'\s+', engine='python', dtype={0: str, 1: str})
    d.columns = [c.lstrip('#') for c in d.columns]
    return d


def main():
    args = parse_args()
    ph = read_table(args.pheno, 'phenotype')
    cv = read_table(args.cov, 'covariates')
    idc = args.id_col

    for name, d in (('phenotype', ph), ('covariates', cv)):
        if idc not in d.columns:
            raise SystemExit(f'{name} file has no {idc} column (has: {list(d.columns)})')
    if args.pheno_col not in ph.columns:
        raise SystemExit(f'phenotype file has no {args.pheno_col} column')

    ph[idc] = ph[idc].astype(str)
    cv[idc] = cv[idc].astype(str)
    # FID is carried by both files and is not a covariate; drop the duplicate.
    cv = cv.drop(columns=[c for c in ('FID',) if c in cv.columns])
    m = ph[[idc, args.pheno_col] + [c for c in ('FID',) if c in ph.columns]].merge(
        cv, on=idc, how='inner')

    n_pheno, n_cov, n_join = len(ph), len(cv), len(m)
    if args.keep_fam:
        fam = pd.read_csv(args.keep_fam, sep=r'\s+', header=None, usecols=[1],
                          names=[idc], dtype={1: str})
        m = m[m[idc].isin(set(fam[idc].astype(str)))]
    n_out = len(m)

    # SAIGE reads a binary phenotype as 0/1. plink's convention is 2 = case,
    # 1 = control, so recode explicitly rather than letting SAIGE guess.
    v = pd.to_numeric(m[args.pheno_col], errors='coerce')
    if set(v.dropna().unique()) <= {1.0, 2.0}:
        m[args.pheno_col] = v.map({1.0: 0, 2.0: 1}).astype('Int64')
    m = m[m[args.pheno_col].notna()]

    m.to_csv(args.out, sep='\t', index=False, na_rep='NA')
    n_case = int((m[args.pheno_col] == 1).sum())
    print(f'[merge_pheno_cov] pheno {n_pheno:,} x cov {n_cov:,} -> join {n_join:,} '
          f'-> in genotypes {n_out:,} ({n_case:,} cases / {n_out - n_case:,} controls), '
          f'{m.shape[1]} columns -> {args.out}')
    if n_out == 0:
        raise SystemExit('no samples left after the join — check the ID conventions')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in merge_pheno_cov: {e}', file=sys.stderr)
        sys.exit(1)
