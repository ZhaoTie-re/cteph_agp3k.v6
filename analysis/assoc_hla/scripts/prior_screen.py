#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Test THIS study against what has already been published for CTEPH in
#           Japanese patients, on the alleles those reports named and in the
#           stratum they named them in.
#
#           WHY THIS IS NOT PART OF THE SCAN. The scan asks "what is associated
#           here" and pays a multiple-testing price for asking it of every
#           marker. This asks a different question -- "does the published claim
#           hold in these data" -- of a handful of pre-specified alleles, with a
#           direction and a magnitude stated in advance. Two or three tests, no
#           discovery, so the threshold is 0.05 and it is not competing with the
#           scan for significance.
#
#           THE STRATUM IS THE POINT. Kominami et al. (J Hum Genet 2009;54:108-114)
#           found HLA-B*52:01 (OR 2.47) and HLA-DPB1*02:02 (OR 5.07) in CTEPH
#           patients WITHOUT deep vein thrombosis, and explicitly not in those
#           with it. Testing the whole case series against controls therefore
#           tests something the paper did not claim, and would dilute a real
#           effect if one existed.
#
#           THIS TEST IS FIXED-EFFECTS, IN BOTH SENSES. It refits a logistic
#           model on a SUBSET of the cases, and the GRM null model the random
#           arm fits is estimated on the whole case series -- it does not apply
#           to a subset. So there is one stratified screen, not one per model,
#           and it says so rather than emitting a duplicate row under a label it
#           has not earned.
#
#           POWER IS REPORTED, ALWAYS. A null result is only informative if the
#           study could have seen the effect. Every row carries the power to
#           detect the REPORTED odds ratio at alpha 0.05, and the smallest odds
#           ratio this stratum could detect with 80 % power, so "did not
#           replicate" can be told apart from "could not have".
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
    p.add_argument('--raw', required=True, help='plink2 --export A over the tested alleles')
    p.add_argument('--pheno-covar', required=True)
    p.add_argument('--subgroups', required=True, help='prep_subgroups.py output')
    p.add_argument('--prior-reports', required=True, help='info/prior_reports.tsv')
    p.add_argument('--marker-qc', required=True)
    p.add_argument('--covars', required=True)
    p.add_argument('--pheno-name', default='PHENO1')
    p.add_argument('--sample-id-col', default='IID')
    p.add_argument('--cohort', required=True)
    p.add_argument('--model', default='fixed_stratified',
                   help='recorded in the output. This test is ALWAYS fixed-effects: it '
                        'refits a logistic model on a case SUBSET, and a GRM null fitted '
                        'on the full case series does not apply to a subset. Emitting a '
                        'row labelled "random" that is not a random-model fit would be a '
                        'quiet falsehood.')
    p.add_argument('--alpha', type=float, default=0.05)
    p.add_argument('--out', default='prior_screen.tsv')
    return p.parse_args()


def power_at(n_case, n_ctrl, af, odds_ratio, alpha=0.05):
    """Power of the additive Wald test to detect `odds_ratio` at this allele
    frequency and these group sizes. A normal approximation, which is what the
    published power calculations for this design use too."""
    if not (0 < af < 1) or not np.isfinite(odds_ratio) or odds_ratio <= 0 or odds_ratio == 1:
        return np.nan
    p1 = af * odds_ratio / (1 + af * (odds_ratio - 1))
    n1, n0 = 2 * n_case, 2 * n_ctrl
    if min(n1 * p1 * (1 - p1), n0 * af * (1 - af)) <= 0:
        return np.nan
    se = np.sqrt(1 / (n1 * p1 * (1 - p1)) + 1 / (n0 * af * (1 - af)))
    z = abs(np.log(odds_ratio)) / se
    zc = stats.norm.isf(alpha / 2)
    return float(stats.norm.sf(zc - z) + stats.norm.cdf(-zc - z))


def min_detectable_or(n_case, n_ctrl, af, alpha=0.05, target=0.80):
    lo, hi = 1.0 + 1e-6, 20.0
    if power_at(n_case, n_ctrl, af, hi, alpha) < target:
        return np.nan
    for _ in range(60):
        mid = (lo + hi) / 2
        if power_at(n_case, n_ctrl, af, mid, alpha) < target:
            lo = mid
        else:
            hi = mid
    return float(hi)


def main():
    args = parse_args()
    import statsmodels.api as sm

    raw = pd.read_csv(args.raw, sep=r'\s+', dtype={'FID': str, 'IID': str}).set_index('IID')
    pc = pd.read_csv(args.pheno_covar, sep='\t',
                     dtype={args.sample_id_col: str}).set_index(args.sample_id_col)
    sub = pd.read_csv(args.subgroups, sep='\t',
                      dtype={args.sample_id_col: str}).set_index(args.sample_id_col)
    pri = pd.read_csv(args.prior_reports, sep='\t')
    qc = pd.read_csv(args.marker_qc, sep='\t')

    covars = [c.strip() for c in args.covars.split(',') if c.strip()]
    ids = [i for i in raw.index if i in pc.index and i in sub.index]
    y = pd.to_numeric(pc.loc[ids, args.pheno_name], errors='coerce')
    X0 = pc.loc[ids, covars].apply(pd.to_numeric, errors='coerce')
    grp = sub.loc[ids, 'dvt_status']

    colmap = {}
    for c in raw.columns:
        mid, _, allele = c.rpartition('_')
        if allele in ('A', 'P') and mid:
            colmap[mid] = (c, allele)

    strata = {'all': y.notna(),
              'dvt_negative': (y == 0) | ((y == 1) & (grp == 'dvt_negative')),
              'dvt_positive': (y == 0) | ((y == 1) & (grp == 'dvt_positive'))}

    rows = []
    for r in pri.itertuples():
        mid = r.marker
        base = {'cohort': args.cohort, 'model': args.model, 'marker': mid,
                'source': r.source, 'reported_stratum': r.reported_stratum,
                'reported_or': r.reported_or}
        if not int(r.testable):
            rows.append({**base, 'stratum': r.reported_stratum, 'status': 'not_testable',
                         'note': r.note})
            continue
        if mid not in colmap:
            q = qc[qc.id == mid]
            why = ('dropped by QC: ' + str(q.drop_reason.iloc[0])) if len(q) \
                else 'no such marker in this cohort'
            rows.append({**base, 'stratum': r.reported_stratum, 'status': why, 'note': r.note})
            continue
        col, allele = colmap[mid]
        g = pd.to_numeric(raw.loc[ids, col], errors='coerce')
        if allele == 'A':
            g = 2 - g

        for sname, mask in strata.items():
            m = mask & y.notna() & g.notna() & X0.notna().all(axis=1)
            n_case = int((y[m] == 1).sum())
            n_ctrl = int((y[m] == 0).sum())
            af = float(g[m].sum() / (2 * m.sum())) if m.sum() else np.nan
            rec = {**base, 'stratum': sname, 'n_case': n_case, 'n_control': n_ctrl,
                   'af': round(af, 6) if np.isfinite(af) else np.nan, 'status': 'ok'}
            try:
                fit = sm.Logit(y[m], sm.add_constant(
                    pd.concat([X0[m], g[m].rename('G')], axis=1))).fit(disp=0, maxiter=200)
                b, se = float(fit.params['G']), float(fit.bse['G'])
                # 'or' is a Python keyword: as a column name it makes
                # itertuples rename it positionally, which breaks silently.
                rec.update({'odds_ratio': round(float(np.exp(b)), 4),
                            'l95': round(float(np.exp(b - 1.959963984540054 * se)), 4),
                            'u95': round(float(np.exp(b + 1.959963984540054 * se)), 4),
                            'p': float(fit.pvalues['G'])})
            except Exception as e:                    # noqa: BLE001
                rec['status'] = f'fit_failed:{type(e).__name__}'
            ror = pd.to_numeric(pd.Series([r.reported_or]), errors='coerce').iloc[0]
            rec['power_for_reported_or'] = (round(power_at(n_case, n_ctrl, af, ror), 4)
                                            if np.isfinite(ror) else np.nan)
            rec['min_or_80pct_power'] = round(min_detectable_or(n_case, n_ctrl, af), 4) \
                if np.isfinite(af) else np.nan
            # The strongest statement a null makes: does our interval EXCLUDE the
            # published point estimate? That is a different claim from "P > 0.05".
            rec['reported_or_outside_ci'] = (
                int(not (rec.get('l95', np.nan) <= ror <= rec.get('u95', np.nan)))
                if np.isfinite(ror) and 'l95' in rec else pd.NA)
            rec['note'] = r.note
            rows.append(rec)

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep='\t', index=False, na_rep='NA')

    print(f'[prior_screen] {args.cohort}/{args.model}: {len(pri)} prior claim(s) '
          f'-> {args.out}')
    for r in out[out.status == 'ok'].itertuples():
        if r.stratum != r.reported_stratum:
            continue
        pw = f'{r.power_for_reported_or:.0%}' if pd.notna(r.power_for_reported_or) else 'n/a'
        excl = (' CI EXCLUDES the reported OR'
                if pd.notna(r.reported_or_outside_ci) and r.reported_or_outside_ci == 1
                else '')
        print(f'    {r.marker:16s} [{r.stratum:12s}] n_case={int(r.n_case):4d}  '
              f'OR={r.odds_ratio:5.2f} [{r.l95:.2f}-{r.u95:.2f}] P={r.p:.3f}  '
              f'reported {r.reported_or}  power {pw}{excl}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in prior_screen: {e}', file=sys.stderr)
        sys.exit(1)
