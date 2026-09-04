#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Clinical subgroups for the case series, from the study's own sample
#           sheet. At present that means DVT status, which is not a covariate but
#           a STRATUM: the only prior Japanese CTEPH HLA report found its
#           associations exclusively in patients WITHOUT deep vein thrombosis
#           (Kominami et al., J Hum Genet 2009), so testing the whole case series
#           against controls tests a different hypothesis from the published one.
#
#           WGS_ID IS THE JOIN KEY. It matches the genotype .fam exactly --
#           3,101 of 3,101 on the full cohort. JointCall_Sample_Code and
#           RADDAR_J_ID match nothing, and CGM_ID matches 63; using any of them
#           would silently drop almost every sample.
#
#           A CONTROL HAS NO DVT STATUS, and that is not missing data. The field
#           is only recorded for the PH series, so controls are marked
#           'not_applicable' rather than left blank, which keeps a stratified
#           comparison from quietly dropping them.
# Component: assoc_hla
# ---------------------------------------------------------------------------
import argparse
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample-sheet', required=True, help='info/*.xlsx')
    p.add_argument('--sheet', default='Sheet1')
    p.add_argument('--pheno-covar', required=True, help='PREP_PHENO_COV output')
    p.add_argument('--id-col', default='WGS_ID')
    p.add_argument('--sample-id-col', default='IID')
    p.add_argument('--pheno-name', default='PHENO1')
    p.add_argument('--cohort', required=True)
    p.add_argument('--out', default='subgroups.tsv')
    return p.parse_args()


def main():
    args = parse_args()
    meta = pd.read_excel(args.sample_sheet, sheet_name=args.sheet, dtype=str)
    if args.id_col not in meta.columns:
        raise SystemExit(f'ABORT: {args.sample_sheet} has no {args.id_col!r} column')
    if 'DVT_Presence' not in meta.columns:
        raise SystemExit(f'ABORT: {args.sample_sheet} has no DVT_Presence column')

    pc = pd.read_csv(args.pheno_covar, sep='\t', dtype={args.sample_id_col: str})
    n_match = len(set(meta[args.id_col].dropna()) & set(pc[args.sample_id_col]))
    if n_match < 0.99 * len(pc):
        raise SystemExit(f'ABORT: {args.id_col} matches only {n_match} of {len(pc)} cohort '
                         f'samples. The join key is wrong -- WGS_ID is the one that matches.')

    dvt = dict(zip(meta[args.id_col], meta.DVT_Presence))
    out = pc[[args.sample_id_col, args.pheno_name]].copy()
    raw = out[args.sample_id_col].map(dvt)
    is_case = pd.to_numeric(out[args.pheno_name], errors='coerce') == 1

    status = pd.Series('not_applicable', index=out.index)
    status[is_case & (raw == 'Yes')] = 'dvt_positive'
    status[is_case & (raw == 'No')] = 'dvt_negative'
    status[is_case & ~raw.isin(['Yes', 'No'])] = 'dvt_unknown'
    out['dvt_status'] = status
    out['cohort'] = args.cohort
    out.to_csv(args.out, sep='\t', index=False)

    n = out.dvt_status.value_counts()
    n_case = int(is_case.sum())
    print(f'[prep_subgroups] {args.cohort}: {len(out):,} sample(s), {n_case:,} case(s) '
          f'-> {args.out}')
    for k in ('dvt_negative', 'dvt_positive', 'dvt_unknown'):
        v = int(n.get(k, 0))
        print(f'    {k:14s} {v:4d}' + (f'  ({v / n_case:.1%} of cases)' if n_case else ''))
    print(f'    not_applicable {int(n.get("not_applicable", 0)):4d}  (controls)')
    if int(n.get('dvt_negative', 0)) < 30:
        print('    WARNING: the DVT-negative stratum is small; a non-replication there '
              'may be a power statement rather than a result.', file=sys.stderr)


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in prep_subgroups: {e}', file=sys.stderr)
        sys.exit(1)
