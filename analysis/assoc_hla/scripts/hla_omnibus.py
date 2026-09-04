#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The amino-acid OMNIBUS test -- the layer that makes this HLA fine
#           mapping rather than a marker scan. It asks, at each position, whether
#           the residue CONTENT of that position associates with disease, rather
#           than whether one particular residue does.
#
#           THE TEST. At a position carrying m residues, the full model adds m-1
#           residue dosages to the covariate-only model (one is dropped: the
#           dosages sum to the position's determined chromosome count, so all m
#           would be collinear). The statistic is the deviance
#               D = 2 * (logL_full - logL_null)
#           compared to chi-squared with m-1 degrees of freedom. This is the test
#           Hirata et al. Nat. Genet. 2019 used on the Japanese MHC, and the one
#           Nat. Protoc. 2023 (doi:10.1038/s41596-023-00853-4) specifies.
#
#           WHY IT IS DONE HERE AND NOT BY THE ENGINES. plink2 and SAIGE both
#           test one marker at a time. A joint m-1 df test over a residue set is
#           not something either exposes, so the fixed-model omnibus is fitted
#           directly. The random-model counterpart is SAIGE-GENE+'s group test,
#           which is a DIFFERENT test (SKAT-O over the residue set, not a
#           likelihood ratio) and is reported beside this one, never merged with
#           it.
#
#           SEPARATION. A residue carried by very few cases drives the logistic
#           fit to the boundary and the deviance to infinity. Fits that do not
#           converge, or whose deviance is not finite, are reported with
#           status != 'ok' and are never silently dropped.
# Component: assoc_hla
# ---------------------------------------------------------------------------
import argparse
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings('ignore')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw', required=True,
                   help='plink2 --export A output over the tested residue markers')
    p.add_argument('--pheno-covar', required=True, help='PREP_PHENO_COV output')
    p.add_argument('--marker-qc', required=True, help='hla_marker_qc.py output')
    p.add_argument('--covars', required=True, help='comma-separated covariate columns')
    p.add_argument('--pheno-name', default='PHENO1')
    p.add_argument('--sample-id-col', default='IID')
    p.add_argument('--cohort', required=True)
    p.add_argument('--model', default='fixed')
    p.add_argument('--condition-on', default=None,
                   help='comma-separated marker ids added to BOTH models, for the '
                        'conditional rounds')
    p.add_argument('--out', default='omnibus.tsv')
    return p.parse_args()


def counted_allele_columns(raw):
    """plink2 writes <ID>_<counted allele>. Our markers are A/P, so the suffix is
    one character and the id is everything before the final underscore. The
    counted allele is recorded so the sign convention is explicit, not assumed."""
    out = {}
    for c in raw.columns:
        if c in ('FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'):
            continue
        mid, _, allele = c.rpartition('_')
        if allele in ('A', 'P') and mid:
            out[mid] = (c, allele)
    return out


def main():
    args = parse_args()

    raw = pd.read_csv(args.raw, sep=r'\s+', dtype={'FID': str, 'IID': str})
    pc = pd.read_csv(args.pheno_covar, sep='\t', dtype={args.sample_id_col: str})
    qc = pd.read_csv(args.marker_qc, sep='\t')

    covars = [c.strip() for c in args.covars.split(',') if c.strip()]
    missing = [c for c in covars + [args.pheno_name] if c not in pc.columns]
    if missing:
        raise SystemExit(f'ABORT: {args.pheno_covar} lacks {missing}')

    colmap = counted_allele_columns(raw)
    dat = raw.set_index('IID')
    pc = pc.set_index(args.sample_id_col)
    ids = [i for i in dat.index if i in pc.index]
    if len(ids) < len(dat):
        print(f'[hla_omnibus] {len(dat) - len(ids)} genotyped sample(s) absent from '
              f'{args.pheno_covar}', file=sys.stderr)
    dat, pc = dat.loc[ids], pc.loc[ids]

    y = pd.to_numeric(pc[args.pheno_name], errors='coerce')
    # plink .fam / plink2 --pheno both use 1 = control, 2 = case; the model_inputs
    # phenotype file uses 0/1. Normalise to 0/1 rather than assuming either.
    vals = set(y.dropna().unique())
    if vals <= {1.0, 2.0}:
        y = y - 1
    elif not vals <= {0.0, 1.0}:
        raise SystemExit(f'ABORT: {args.pheno_name} takes values {sorted(vals)}; expected '
                         f'0/1 or 1/2')

    X0 = pc[covars].apply(pd.to_numeric, errors='coerce')
    cond_ids = [c.strip() for c in (args.condition_on or '').split(',') if c.strip()]
    for cid in cond_ids:
        if cid not in colmap:
            raise SystemExit(f'ABORT: cannot condition on {cid!r}: not in {args.raw}')
        col, allele = colmap[cid]
        v = dat[col]
        X0[f'COND_{cid}'] = (2 - v) if allele == 'A' else v

    import statsmodels.api as sm

    elig = qc[qc.omnibus_eligible.astype(str).str.lower().isin(('true', '1'))]
    rows = []
    for pos, grp in elig.groupby('position'):
        ids_here = [i for i in grp.id if i in colmap]
        m = len(ids_here)
        if m < 2:
            continue
        Z = pd.DataFrame(index=dat.index)
        for mid in ids_here:
            col, allele = colmap[mid]
            v = pd.to_numeric(dat[col], errors='coerce')
            Z[mid] = (2 - v) if allele == 'A' else v
        # Drop ONE residue as the reference level. The dosages at a position sum
        # to the number of determined chromosomes, so the full set is collinear
        # with the intercept. The commonest residue is the reference, which is
        # the conventional choice and makes the remaining effects readable.
        ref_residue = Z.sum().idxmax()
        Zr = Z.drop(columns=[ref_residue])

        ok = y.notna() & X0.notna().all(axis=1) & Z.notna().all(axis=1)
        n = int(ok.sum())
        rec = {'cohort': args.cohort, 'model': args.model, 'position': pos,
               'gene': grp.gene.iloc[0], 'n_residues': m, 'df': m - 1, 'n': n,
               'reference_residue': ref_residue,
               'residues': ';'.join(ids_here),
               'determined_rate': grp.determined_rate.iloc[0],
               'deviance': np.nan, 'P': np.nan, 'status': 'ok'}
        if n < 50 or Zr.loc[ok].std().min() == 0:
            rec['status'] = 'too_few' if n < 50 else 'monomorphic_in_complete_cases'
            rows.append(rec)
            continue
        try:
            base = sm.add_constant(X0.loc[ok], has_constant='add')
            f0 = sm.Logit(y[ok], base).fit(disp=0, maxiter=200)
            f1 = sm.Logit(y[ok], pd.concat([base, Zr.loc[ok]], axis=1)).fit(disp=0, maxiter=200)
            dev = 2.0 * (f1.llf - f0.llf)
            if not (f0.mle_retvals.get('converged') and f1.mle_retvals.get('converged')):
                rec['status'] = 'not_converged'
            elif not np.isfinite(dev) or dev < 0:
                rec['status'] = 'bad_deviance'
            else:
                rec['deviance'] = round(float(dev), 6)
                rec['P'] = float(stats.chi2.sf(dev, m - 1))
        except Exception as e:                       # noqa: BLE001
            rec['status'] = f'error:{type(e).__name__}'
        rows.append(rec)

    out = pd.DataFrame(rows).sort_values('P', na_position='last')
    out.to_csv(args.out, sep='\t', index=False, na_rep='NA')

    good = out[out.status == 'ok']
    print(f'[hla_omnibus] {args.cohort}/{args.model}: {len(out):,} position(s) tested, '
          f'{len(good):,} with a usable fit -> {args.out}')
    if cond_ids:
        print(f'    conditioned on: {", ".join(cond_ids)}')
    if len(good):
        b = good.iloc[0]
        print(f'    strongest: {b.position} ({b.gene}) df={b.df} '
              f'deviance={b.deviance:.2f} P={b.P:.3e}')
    bad = out[out.status != 'ok'].status.value_counts()
    if len(bad):
        print('    not usable: ' + ', '.join(f'{k} {v}' for k, v in bad.items()))


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in hla_omnibus: {e}', file=sys.stderr)
        sys.exit(1)
