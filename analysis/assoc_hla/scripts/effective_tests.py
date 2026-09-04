#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : How many INDEPENDENT tests this scan really performed, and therefore
#           what "significant" means here.
#
#           WHY NOT JUST 5e-8. That threshold is calibrated for a genome-wide
#           scan of a million-odd independent common variants. This scan tests
#           about 150 alleles and 1,000 residues inside a 3.4 Mb region whose LD
#           is the most extreme in the genome -- the residues at one position are
#           nearly collinear by construction, and whole haplotypes travel
#           together across genes. Applying 5e-8 here is not conservative, it is
#           mis-calibrated: it answers a question about a different experiment.
#
#           WHY NOT PLAIN BONFERRONI EITHER. 0.05/1000 assumes 1,000 independent
#           tests. They are nothing of the kind.
#
#           WHAT IS DONE. Li & Ji (2005, Heredity 95:221-227) estimate the
#           effective number of independent tests from the eigenvalues of the
#           marker correlation matrix:
#               Meff = sum_i [ I(lambda_i >= 1) + (lambda_i - floor(lambda_i)) ]
#           and the threshold is 0.05 / Meff. All three numbers -- Bonferroni,
#           effective, and the genome-wide convention -- are written out, with
#           the effective one marked as the primary call, so a reader can apply
#           whichever their journal expects without recomputing anything.
# Component: assoc_hla
# ---------------------------------------------------------------------------
import argparse
import sys

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode', default='marker', choices=('marker', 'omnibus'),
                   help="'marker' scores the single-marker scan; 'omnibus' scores the "
                        "position-level amino-acid tests, which are FEWER and differently "
                        "correlated and therefore need their own threshold")
    p.add_argument('--raw', required=True, help='plink2 --export A over the tested markers')
    p.add_argument('--sumstat', help='marker mode: the sumstats table')
    p.add_argument('--omnibus', help='omnibus mode: hla_omnibus.py output')
    p.add_argument('--marker-qc', help='omnibus mode: to map residues to their position')
    p.add_argument('--cohort', required=True)
    p.add_argument('--model', required=True)
    p.add_argument('--marker-class', required=True)
    p.add_argument('--alpha', type=float, default=0.05)
    p.add_argument('--p-genomewide', type=float, default=5e-8)
    p.add_argument('--out', default='significance.tsv')
    return p.parse_args()


def li_ji_meff(corr):
    """Li & Ji (2005). The integer part counts an eigenvalue >= 1 as one full
    test; the fractional part credits partial independence."""
    ev = np.linalg.eigvalsh(corr)
    ev = np.clip(ev, 0.0, None)
    return float(np.sum((ev >= 1.0).astype(float) + (ev - np.floor(ev))))


def lambda_gc(p):
    from scipy import stats
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    if not len(p):
        return float('nan')
    return float(stats.chi2.isf(np.median(p), df=1) / 0.4549364231195729)


def position_matrix(G, qc):
    """One column per amino-acid position: the first principal component of that
    position's residue dosages.

    A position-level test has no single genotype vector to correlate -- it is a
    joint test over m residues -- so the effective-test count needs a stand-in for
    "this position's genotype". PC1 of the position's own residue dosages is that:
    it is the dominant contrast the omnibus is testing, and correlations between
    positions' PC1s carry the LD between positions. This is an APPROXIMATION and
    is labelled as one; the Bonferroni count beside it needs no such assumption.
    """
    cols_by_pos = {}
    for _i, r in qc.iterrows():
        if r.get('marker_class') != 'residue':
            continue
        for c in G.columns:
            if c.rpartition('_')[0] == r['id']:
                cols_by_pos.setdefault(r['position'], []).append(c)
    out = {}
    for pos, cols in cols_by_pos.items():
        if len(cols) < 2:
            continue
        X = G[cols].to_numpy(dtype=float)
        X = np.where(np.isnan(X), np.nanmean(X, axis=0), X)
        X = X - X.mean(axis=0)
        sd = X.std(axis=0)
        if not np.all(np.isfinite(sd)) or sd.max() <= 0:
            continue
        try:
            _u, _s, vt = np.linalg.svd(X, full_matrices=False)
            out[pos] = X @ vt[0]
        except np.linalg.LinAlgError:
            continue
    return pd.DataFrame(out)


def main():
    args = parse_args()
    raw = pd.read_csv(args.raw, sep=r'\s+', dtype={'FID': str, 'IID': str})
    drop = {'FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'}
    G = raw[[c for c in raw.columns if c not in drop]].apply(pd.to_numeric, errors='coerce')

    if args.mode == 'omnibus':
        if not (args.omnibus and args.marker_qc):
            raise SystemExit('ABORT: --mode omnibus needs --omnibus and --marker-qc')
        om = pd.read_csv(args.omnibus, sep='\t')
        usable = om[om.status == 'ok'].rename(columns={'P': 'P'})
        qc = pd.read_csv(args.marker_qc, sep='\t')
        qc = qc[qc.tested.astype(str).str.lower().isin(('true', '1'))]
        G = position_matrix(G, qc)
        G = G[[c for c in G.columns if c in set(usable.position)]]
    else:
        if not args.sumstat:
            raise SystemExit('ABORT: --mode marker needs --sumstat')
        ss = pd.read_csv(args.sumstat, sep='\t')
        usable = ss[ss.ERRCODE == '.']
    # Markers that do not vary carry no information and make the correlation
    # matrix singular; they are counted as zero independent tests, not dropped
    # silently from the denominator.
    sd = G.std()
    G = G[[c for c in G.columns if sd.get(c, 0) > 0]]
    n_const = int((sd <= 0).sum())

    corr = np.corrcoef(G.fillna(G.mean()).to_numpy().T)
    corr = np.nan_to_num(corr, nan=0.0)
    np.fill_diagonal(corr, 1.0)
    meff = li_ji_meff(corr)

    n_tested = len(usable)
    rows = [
        {'threshold': 'effective_tests', 'value': args.alpha / meff if meff > 0 else np.nan,
         'basis': (f'Li & Ji Meff = {meff:.2f} of {G.shape[1]} '
                   + ('positions, via each position\'s PC1' if args.mode == 'omnibus'
                      else 'varying markers')),
         'primary': 1},
        {'threshold': 'bonferroni', 'value': args.alpha / n_tested if n_tested else np.nan,
         'basis': (f'{args.alpha} / {n_tested} '
                   + ('positions with a usable fit' if args.mode == 'omnibus'
                      else 'tested markers')), 'primary': 0},
        {'threshold': 'genome_wide', 'value': args.p_genomewide,
         'basis': 'the genome-wide convention, for comparability with the SNP scans',
         'primary': 0},
    ]
    out = pd.DataFrame(rows)
    out['cohort'], out['model'] = args.cohort, args.model
    out['marker_class'] = args.marker_class
    out['mode'] = args.mode
    out['n_tested'] = n_tested
    out['n_markers_constant'] = n_const
    out['meff'] = round(meff, 4)
    # lambda_GC assumes every test has 1 df. The omnibus does not -- its tests
    # carry m-1 df and m varies by position -- so the statistic is not defined
    # there and is left blank rather than printed as a number that looks like a
    # calibration check and is not one.
    out['lambda_gc'] = (round(lambda_gc(usable.P), 4) if args.mode == 'marker'
                        else pd.NA)
    for r in rows:
        thr = r['value']
        out.loc[out.threshold == r['threshold'], 'n_significant'] = \
            int((usable.P < thr).sum()) if np.isfinite(thr) else 0
    out['n_significant'] = out.n_significant.astype(int)
    out.to_csv(args.out, sep='\t', index=False, na_rep='NA')

    print(f'[effective_tests] {args.cohort}/{args.model}/{args.marker_class}: '
          f'{n_tested:,} tested, Meff = {meff:.2f} '
          f'({meff / max(G.shape[1], 1):.1%} of {G.shape[1]:,} varying markers)')
    if args.mode == 'marker':
        print(f'    lambda_GC = {out.lambda_gc.iloc[0]:.3f}')
    else:
        print('    lambda_GC not defined: the omnibus tests carry m-1 df, not 1')
    for r in out.itertuples():
        star = ' <- primary' if r.primary else ''
        print(f'    {r.threshold:16s} P < {r.value:.3e}  ->  {r.n_significant:3d} '
              f'significant{star}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in effective_tests: {e}', file=sys.stderr)
        sys.exit(1)
