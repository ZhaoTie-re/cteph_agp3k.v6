#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The phenotype and covariates, on exactly the genotyped samples, in
#           the coding rvtest needs, from ONE normalised vector -- plus the
#           case/control keep lists the per-gene counts are taken over.
#
#           PHENOTYPE CODING IS NORMALISED, NOT ASSUMED -- and here that is not
#           defensive style, it is the fix for a live defect:
#
#             * model_inputs/<cohort>/phenotype/pheno.tsv is 0/1 today. It was
#               1/2 when the earlier gene scan ran and was rewritten in place
#               later; the manifest beside it still claims 1/2.
#             * rvtests/src/DataLoader.cpp:1248 reads 0 as MISSING, not control.
#               A 0/1 column still passes _isBinaryPhenotype (values all in
#               {0,1,2}), maps 0 -> -1, then drops every negative row.
#
#           Composed, those give: hundreds of samples, all controls, zero cases,
#           `WARN There are no case!`, and EXIT 0. Every downstream P is then
#           garbage that looks fine -- plausible magnitudes, plausible lambda, a
#           plausible Manhattan. Nothing downstream can tell it from a null scan,
#           which is why this is fixed here and gated inside RVTEST.
#
#           rvtest wants 1/2. pheno_coding.json carries the counts the RVTEST
#           post-flight gate checks the engine's own log against, and the two
#           keep lists are written from the same vector so the per-gene allele
#           and carrier counts can never disagree with the association about
#           who is a case.
# Component: assoc_gene
# ---------------------------------------------------------------------------
import argparse
import json
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pheno', required=True, help='model_inputs phenotype/pheno.tsv')
    p.add_argument('--covar', required=True, help='model_inputs covariates/<pc>.sex.tsv')
    p.add_argument('--fam', required=True,
                   help='the tuning.rv filtered .fam -- DEFINES the sample set and its order')
    p.add_argument('--cohort', required=True,
                   help='cohort tag, written into pheno_coding.json so COLLECT_ALL can key '
                        'the sample counts without parsing a staged path')
    p.add_argument('--pheno-col', default='PHENO1')
    p.add_argument('--id-col', default='IID')
    p.add_argument('--out-rvt-pheno', default='pheno_rvt.tsv')
    p.add_argument('--out-rvt-covar', default='covar_rvt.tsv')
    p.add_argument('--out-case-keep', default='case.keep')
    p.add_argument('--out-control-keep', default='control.keep')
    p.add_argument('--out-json', default='pheno_coding.json')
    return p.parse_args()


def read_tsv(path, id_col):
    d = pd.read_csv(path, sep=r'\s+', dtype=str)
    d.columns = [c.lstrip('#') for c in d.columns]
    if id_col not in d.columns:
        raise SystemExit(f'ABORT: {path} has no {id_col!r} column; found {list(d.columns)[:6]}')
    return d


def normalise(y, col):
    """Accept 0/1 or 1/2; return 0/1 (0 = control, 1 = case) and the mapping used."""
    vals = set(y.dropna().unique())
    if vals <= {1.0, 2.0}:
        return (y - 1).astype('Int64'), '1/2 -> 0/1'
    if vals <= {0.0, 1.0}:
        return y.astype('Int64'), '0/1, unchanged'
    raise SystemExit(f'ABORT: {col} takes {sorted(vals)}; expected 0/1 or 1/2. '
                     f'This component refuses to guess which class is which.')


def main():
    args = parse_args()
    ph = read_tsv(args.pheno, args.id_col)
    cv = read_tsv(args.covar, args.id_col)
    fam = pd.read_csv(args.fam, sep=r'\s+', header=None, dtype=str,
                      names=['FID', 'IID', 'PAT', 'MAT', 'SEX_FAM', 'PHENO_FAM'])

    if args.pheno_col not in ph.columns:
        raise SystemExit(f'ABORT: {args.pheno} has no {args.pheno_col!r}')

    # The .fam defines the sample set AND the order.
    d = fam[['FID', 'IID']].merge(ph[[args.id_col, args.pheno_col]],
                                  left_on='IID', right_on=args.id_col, how='left')
    if args.id_col != 'IID' and args.id_col in d.columns:
        d = d.drop(columns=[args.id_col])
    cv2 = cv.drop(columns=[c for c in ('FID',) if c in cv.columns])
    d = d.merge(cv2, left_on='IID', right_on=args.id_col, how='left')
    if args.id_col != 'IID' and args.id_col in d.columns:
        d = d.drop(columns=[args.id_col])

    y01, coding = normalise(pd.to_numeric(d[args.pheno_col], errors='coerce'), args.pheno_col)
    d[args.pheno_col] = y01

    n_missing = int(d[args.pheno_col].isna().sum())
    if n_missing:
        raise SystemExit(
            f'ABORT: {n_missing} genotyped sample(s) have no phenotype. The cohorts are '
            f'DEFINED from these files, so this is a join failure -- most likely an id '
            f'type coercion that dropped leading zeros -- not a real gap.')

    n_case = int((d[args.pheno_col] == 1).sum())
    n_ctrl = int((d[args.pheno_col] == 0).sum())
    if n_case == 0 or n_ctrl == 0:
        raise SystemExit(
            f'ABORT: {n_case} case / {n_ctrl} control after normalisation ({coding}). '
            f'An empty class is exactly what the 0/1-into-rvtest defect produces, and '
            f'rvtest would run to completion on it.')

    if 'SEX' not in d.columns:
        raise SystemExit(f'ABORT: {args.covar} has no SEX column; the model uses it.')

    front = ['FID', 'IID', args.pheno_col]
    rest = [c for c in d.columns if c not in front]

    # ---- rvtest: 1/2, lower-case headers, tab -----------------------------
    # PED-like column order; rvtest reads --pheno-name against the lower-cased
    # header, and --covar-name likewise.
    rvp = pd.DataFrame({
        'fid': d['FID'], 'iid': d['IID'],
        'fatid': '0', 'matid': '0',
        'sex': d['SEX'].astype(str),
        args.pheno_col.lower(): (d[args.pheno_col] + 1).astype('Int64'),
    })
    rvp.to_csv(args.out_rvt_pheno, sep='\t', index=False)

    rvc = d[['FID', 'IID'] + [c for c in rest]].copy()
    rvc.columns = [c.lower() for c in rvc.columns]
    rvc.to_csv(args.out_rvt_covar, sep='\t', index=False)

    # ---- the keep lists: FID IID, no header, .fam order --------------------
    d[d[args.pheno_col] == 1][['FID', 'IID']].to_csv(
        args.out_case_keep, sep='\t', index=False, header=False)
    d[d[args.pheno_col] == 0][['FID', 'IID']].to_csv(
        args.out_control_keep, sep='\t', index=False, header=False)

    # ---- the contract the RVTEST post-flight gate checks against ----------
    with open(args.out_json, 'w') as fh:
        json.dump({
            'cohort': args.cohort,
            'n_samples': int(len(d)),
            'n_case': n_case,
            'n_control': n_ctrl,
            'input_coding': coding,
            'rvtest_coding': '1/2 (1=control, 2=case)',
            'pheno_col_rvtest': args.pheno_col.lower(),
            'source_pheno': args.pheno,
            'source_covar': args.covar,
            'source_fam': args.fam,
        }, fh, indent=2)

    print(f'[prep_pheno_cov] {len(d):,} sample(s): {n_case:,} case / {n_ctrl:,} control '
          f'({coding})')
    print(f'    rvtest 1/2 -> {args.out_rvt_pheno} + {args.out_rvt_covar}')
    print(f'    keep lists -> {args.out_case_keep} ({n_case:,}), '
          f'{args.out_control_keep} ({n_ctrl:,})')
    print(f'    contract   -> {args.out_json}')
    print(f'    {len(rest)} covariate column(s) available; the model selects from them')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in prep_pheno_cov: {e}', file=sys.stderr)
        sys.exit(1)
