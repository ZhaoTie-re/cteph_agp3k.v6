#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The numbers of the rare-variant QC table (sheet Rare-variant QC) for
#           one cohort, recomputed from the tuning.rv outputs and checked against
#           them (the computation of the project's rare-variant QC summary):
#             input        samples (case / control) and variants entering the QC
#             sample QC    replay of the outlier rule: Theil-Sen line of the burden
#                          rate on mean depth, robust Z of the residual, |Z| > the
#                          run's threshold; must reproduce the removal list
#             variant QC   delta = |VMISS_30X - VMISS_15X| on the sample-cleaned set;
#                          the count above the knee must reproduce the removal count;
#                          variants left monomorphic by the sample removal, and how
#                          many of them the delta filter removed anyway
#             allele count the minAC sweep and recommendation; the case-control
#                          burden gap at minAC 1 and at the chosen minAC under the
#                          decision model (r ~ case + z(depth); asserted equal to the
#                          sweep) and under the association covariates (r ~ case +
#                          sex + the PCs)
# Inputs  : the tuning.rv files of the cohort (see --help), the sweep, the cohort's
#           sex + PC covariate table
# Outputs : rv_qc.<cohort>.json
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import theilslopes

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gs_common import log  # noqa: E402


def kv(path):
    t = pd.read_csv(path, sep='\t', dtype=str)
    return dict(zip(t.iloc[:, 0], t.iloc[:, 1]))


def gap(d, terms):
    """OLS of the rate on case + terms: the case coefficient, its 95% CI and P, also as % of the mean rate."""
    f = smf.ols('r ~ case + ' + ' + '.join(terms), d).fit()
    lo, hi = f.conf_int().loc['case']
    m = float(d.r.mean())
    b = float(f.params['case'])
    return {'beta': b, 'lo': float(lo), 'hi': float(hi), 'p': float(f.pvalues['case']), 'mean_r': m,
            'pct': 100 * b / m, 'pct_lo': 100 * float(lo) / m, 'pct_hi': 100 * float(hi) / m}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohort', required=True)
    ap.add_argument('--sample-metrics', required=True, help='00.precheck/metrics/sample_metrics.txt.gz')
    ap.add_argument('--outlier-log', required=True, help='00.precheck/01_sample_qc/detect_outlier_samples.log')
    ap.add_argument('--outliers', required=True, help='00.precheck/01_sample_qc/sample_outliers.remove.txt')
    ap.add_argument('--variant-metrics', required=True, help='00.precheck/02_variant_qc/variant_metrics.txt.gz (after the sample QC)')
    ap.add_argument('--depthdiff-stats', required=True, help='00.precheck/02_variant_qc/depthdiff_stats.tsv')
    ap.add_argument('--fam', required=True, help='02.callset_filter/filtered.fam')
    ap.add_argument('--recommendation', required=True, help='05.qc_collect/minac_recommendation.tsv')
    ap.add_argument('--sweep', required=True, help='_comparison/all_cohorts_qc_stats.tsv')
    ap.add_argument('--minac-metrics', nargs='+', required=True, help='<k>=03.qc_metrics/minac<k>/sample_metrics.txt.gz')
    ap.add_argument('--covariates', required=True, help='sex + PC table of the cohort (IID, SEX, PC1_AVG ...)')
    ap.add_argument('--n-pcs', type=int, required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    c = a.cohort
    out = {'cohort': c}

    # ── input and sample QC ─────────────────────────────────────────────────
    m = re.search(r'z_threshold\s*:\s*([\d.]+)', Path(a.outlier_log).read_text())
    assert m, 'z_threshold not in the outlier log'
    zthr = float(m.group(1))
    s = pd.read_csv(a.sample_metrics, sep='\t', dtype={'SampleID': str})
    s['N'] = s.SNumHomRef + s.SNumHet + s.SNumHomAlt
    s['rate'] = 1000 * s.SMinAC / s.N
    out['n_in'] = {'all': int(len(s)), 'case': int((s.Group == 'Case').sum()),
                   'control': int((s.Group == 'Control').sum())}
    slope, intercept = theilslopes(s.rate, s.MeanDP)[:2]
    s['resid'] = s.rate - (slope * s.MeanDP + intercept)
    med, mad = s.resid.median(), (s.resid - s.resid.median()).abs().median()
    s['z'] = 0.6745 * (s.resid - med) / mad
    out_ids = set(s.SampleID[s.z.abs() > zthr])
    removed = set(pd.read_csv(a.outliers, sep='\t', header=None, dtype=str)[1])
    assert out_ids == removed, 'the replay must reproduce the removal list'
    fam = pd.read_csv(a.fam, sep=r'\s+', header=None, dtype={1: str})
    out['n_kept'] = {'all': int(len(fam)), 'case': int((fam[5] == 2).sum()), 'control': int((fam[5] == 1).sum())}
    assert out['n_kept']['all'] == out['n_in']['all'] - len(removed)
    assert not (set(fam[1]) & removed)
    out['sample_qc'] = {'z_threshold': zthr, 'n_removed': len(removed),
                        'removed_case': int(s[s.SampleID.isin(removed)].Group.eq('Case').sum()),
                        'max_abs_z_kept': float(s.z.abs()[~s.SampleID.isin(removed)].max())}

    # ── variant QC on the sample-cleaned set ────────────────────────────────
    dd = kv(a.depthdiff_stats)
    knee = float(dd['knee'])
    n_valid = n_above = n_mono = n_mono_above = 0
    for ch in pd.read_csv(a.variant_metrics, sep='\t', chunksize=2_000_000,
                          usecols=['MinAC', 'VMISS_15X', 'VMISS_30X'],
                          dtype={'VMISS_15X': 'float32', 'VMISS_30X': 'float32'}):
        d = (ch.VMISS_30X - ch.VMISS_15X).abs().to_numpy(np.float64)
        ok = np.isfinite(d)
        n_valid += int(ok.sum())
        n_above += int((d[ok] > knee).sum())
        mono = (ch.MinAC == 0).to_numpy()
        n_mono += int(mono.sum())
        n_mono_above += int((mono & ok & (d > knee)).sum())
    assert n_valid == int(dd['n_valid']) == int(dd['n_total']), (n_valid, dd['n_valid'])
    assert n_above == int(dd['n_removed']) + int(dd['n_protected_rescued']), 'delta > knee must reproduce the removal'
    out['variant_qc'] = {'knee': knee, 'n_in': n_valid, 'n_removed': int(dd['n_removed']),
                         'n_protected_rescued': int(dd['n_protected_rescued']),
                         'n_mono_after_sample_qc': n_mono, 'n_mono_removed_by_delta': n_mono_above}

    # ── the minAC sweep and recommendation ──────────────────────────────────
    rec = kv(a.recommendation)
    out['minac'] = {k: rec[k] for k in ('recommended_minac', 'calibrated', 'nonsig_band', 'beta_group_at_rec',
                                        'p_group_at_rec', 'variant_count_at_rec', 'rule', 'alpha')}
    kappa = int(rec['recommended_minac'])
    sw = pd.read_csv(a.sweep, sep='\t')
    sw = sw[sw.Cohort == c].set_index('MinAC_Threshold')
    vq = out['variant_qc']
    k0, k1, kk = (int(sw.loc[k, 'Variant_Count']) for k in (0, 1, kappa))
    assert k0 == vq['n_in'] - vq['n_removed'], 'the sweep at minAC 0 is the variant-QC output'
    assert k1 == vq['n_in'] - vq['n_removed'] - (vq['n_mono_after_sample_qc'] - vq['n_mono_removed_by_delta']), \
        '... and at minAC 1 that output without the sites the sample removal left monomorphic'
    assert kk == int(rec['variant_count_at_rec'])
    out['sweep'] = {'0': k0, '1': k1, str(kappa): kk}

    # ── the case-control burden gap at minAC 1 and at the chosen minAC ───────
    cov = pd.read_csv(a.covariates, sep='\t', dtype={'IID': str}).set_index('IID')
    pcs = [f'PC{i}_AVG' for i in range(1, a.n_pcs + 1)]
    assert set(['SEX'] + pcs) <= set(cov.columns), cov.columns.tolist()[:15]
    mets = {int(x.split('=', 1)[0]): x.split('=', 1)[1] for x in a.minac_metrics}
    assert {1, kappa} <= set(mets), sorted(mets)
    out['model_check'] = {}
    for k in sorted({1, kappa}):
        d = pd.read_csv(mets[k], sep='\t', dtype={'SampleID': str}).set_index('SampleID')
        d['r'] = 1000 * d.SMinAC / (d.SNumHomRef + d.SNumHet + d.SNumHomAlt)
        d['case'] = (d.Group == 'Case').astype(int)
        d['zdp'] = (d.MeanDP - d.MeanDP.mean()) / d.MeanDP.std()
        n = len(d)
        assert n == out['n_kept']['all'], (k, n)
        d = d.join(cov[['SEX'] + pcs], how='inner')
        assert len(d) == n, 'every sample needs sex and PCs'
        dec, asc = gap(d, ['zdp']), gap(d, ['SEX'] + pcs)
        ref = sw.loc[k]
        assert abs(dec['beta'] - ref.Beta_group_noPC) < 1e-6 and \
            abs(dec['p'] - ref.P_group_noPC) < 1e-4 * max(ref.P_group_noPC, 1e-12) + 1e-9, \
            'the replayed decision model must reproduce the sweep'
        out['model_check'][str(k)] = {'decision': dec, 'association': asc}
    json.dump(out, open(a.out, 'w'), indent=1)
    log(f'{c}: samples {out["n_in"]} -> {out["n_kept"]} (|robust Z| > {zthr:g}: {len(removed)}); variants '
        f'{n_valid:,} -> delta > {knee:g}: {vq["n_removed"]:,}; monomorphic after sample QC {n_mono:,} '
        f'({n_mono_above:,} also above the knee); minAC {kappa}: {kk:,}')


if __name__ == '__main__':
    main()
