#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One engine's output for one (cohort, class) -> the component's one
#           result schema, the family-wise calls, and a one-row scan summary.
#
#           SCHEMA (both models): the marker and its counts in the samples that
#           model analysed (marker_info.tsv, sample_set = model), then
#           obs, or, l95, u95, beta, se, p, p_spa, errcode, test_id, mirror_of,
#           sig_gec, sig_maxt, sig_bonf.
#
#             fixed  : plink2 --glm logistic, no-firth. beta = log(OR); the CI is
#                      plink2's Wald interval; p is the Wald P; errcode is
#                      plink2's ERRCODE ('.' = none). p_spa is empty.
#             random : SAIGE step 2 score test on the reused null model.
#                      beta = Tstat / var (SAIGE's one-step estimate);
#                      se = 1 / sqrt(var), the standard error that matches the
#                      unadjusted score P; or and its 95% CI are exp(beta) and
#                      exp(beta +/- 1.96 se). p = p.value.NA (no saddle-point
#                      approximation, as in the genome-wide SAIGE scan);
#                      p_spa = p.value, kept as a secondary column.
#
#           FAMILY-WISE CALLS, per (case set, cohort, class). The family is the
#           class's tests in the case set (effective_tests.py: polymorphic
#           markers, a two-residue mirror pair counted once), fixed BEFORE any P
#           exists. The three thresholds come from that table and are the same
#           for both models:
#             sig_gec  : p < 0.05 / Me (GEC v0.2)          -- primary
#             sig_maxt : p < the exact max-T threshold      -- validation
#             sig_bonf : p < 0.05 / tests                   -- sensitivity
#           A marker with no P is significant under none. Benjamini-Hochberg is
#           not used: in this region the tests are strongly correlated, so a
#           BH call depends on how many correlated markers move together.
#           lambda_GC is reported for each scan and never applied.
#
#           GATES. The run stops here if the engine analysed other samples than
#           the phenotype file says (the 0/1 coding trap), if a marker is
#           duplicated, if a marker with a minor allele in the analysed samples
#           has no result, or if the effect allele is not 'P'. A marker with no
#           minor allele there (case-only set only) keeps its row with
#           errcode MONOMORPHIC and no P (and is not a test if it is monomorphic
#           in the random-model samples too).
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import json
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

Z975 = stats.norm.ppf(0.975)
CHI2_MEDIAN = stats.chi2.ppf(0.5, 1)
V_THRESHOLDS = ('gec', 'maxt', 'bonf')        # vocab.THRESHOLDS; COLLECT stages no vocab.py


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case-set', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--class', dest='cls', required=True)
    p.add_argument('--model', required=True, choices=('fixed', 'random'))
    p.add_argument('--engine-out', required=True)
    p.add_argument('--engine-log', required=True)
    p.add_argument('--marker-info', required=True)
    p.add_argument('--coding', required=True, help='pheno_coding.json')
    p.add_argument('--n-covariates', type=int, required=True)
    p.add_argument('--thresholds', required=True, help='effective_tests.<class>.tsv')
    p.add_argument('--tests', required=True, help='tests.<class>.tsv')
    p.add_argument('--out-assoc', required=True)
    p.add_argument('--out-qc', required=True)
    return p.parse_args()


def lambda_gc(p):
    p = p[~np.isnan(p)]
    if not len(p):
        return np.nan
    return float(np.median(stats.chi2.isf(p, 1)) / CHI2_MEDIAN)


def read_fixed(path, log):
    g = pd.read_csv(path, sep='\t', dtype={'ID': str})
    if set(g['A1']) != {'P'}:
        sys.exit(f'ERROR: plink2 effect allele is not P for every marker: {sorted(set(g.A1))}')
    out = pd.DataFrame({
        'marker': g['ID'], 'obs': g['OBS_CT'],
        'or': g['OR'], 'l95': g['L95'], 'u95': g['U95'],
        'beta': np.log(g['OR']), 'se': g['LOG(OR)_SE'],
        'p': g['P'], 'p_spa': np.nan, 'errcode': g['ERRCODE'].astype(str)})
    txt = open(log).read()
    m = re.search(r'(\d+) cases and (\d+) controls remaining after main filters', txt)
    c = re.search(r'(\d+) covariates? loaded', txt)
    eng = {'engine_cases': int(m.group(1)) if m else None,
           'engine_controls': int(m.group(2)) if m else None,
           'engine_covariates': int(c.group(1)) if c else None}
    return out, eng


def read_random(path, log):
    s = pd.read_csv(path, sep=r'\s+', dtype={'MarkerID': str})
    if set(s['Allele2']) != {'P'}:
        sys.exit(f'ERROR: SAIGE Allele2 is not P for every marker: {sorted(set(s.Allele2))}')
    var = s['var'].astype(float)
    beta = s['BETA'].astype(float)
    se = 1.0 / np.sqrt(var)
    out = pd.DataFrame({
        'marker': s['MarkerID'], 'obs': s['N_case'] + s['N_ctrl'],
        'or': np.exp(beta), 'l95': np.exp(beta - Z975 * se), 'u95': np.exp(beta + Z975 * se),
        'beta': beta, 'se': se,
        'p': s['p.value.NA'].astype(float), 'p_spa': s['p.value'].astype(float),
        'errcode': np.where(s['p.value.NA'].isna(), 'NO_RESULT', '.')})
    txt = open(log).read()
    nc, nn = s['N_case'].unique(), s['N_ctrl'].unique()
    imp = re.search(r'\$impute_method\s*\n\[1\] "(\w+)"', txt)
    vr = re.search(r'variance Ratio null is\s+([0-9.eE+-]+)', txt)
    eng = {'engine_cases': int(nc[0]) if len(nc) == 1 else None,
           'engine_controls': int(nn[0]) if len(nn) == 1 else None,
           'engine_covariates': None,              # taken from the null model upstream
           'impute_method': imp.group(1) if imp else None,
           'loco_chr6': bool(re.search(r'Leave chromosome\s+6\s+out will be applied', txt)),
           'variance_ratio': float(vr.group(1)) if vr else None,
           'n_spa': int((s['Is.SPA'].astype(str).str.upper() == 'TRUE').sum())}
    return out, eng


def main():
    a = parse_args()
    info = pd.read_csv(a.marker_info, sep='\t', dtype={'marker': str, 'chrom': str})
    info = info[info.sample_set == a.model].drop(columns='sample_set')
    if set(info['case_set']) != {a.case_set}:
        sys.exit(f'ERROR: marker_info is for case set(s) {sorted(set(info.case_set))}, not {a.case_set}')
    coding = json.load(open(a.coding))[a.model]
    m = len(info)
    th = pd.read_csv(a.thresholds, sep='\t').iloc[0]
    tests = pd.read_csv(a.tests, sep='\t', dtype=str, keep_default_na=False)
    if (th['case_set'], th['cohort'], th['class']) != (a.case_set, a.cohort, a.cls):
        sys.exit(f"ERROR: thresholds are for {th['case_set']} {th['cohort']} {th['class']}")

    res, eng = (read_fixed if a.model == 'fixed' else read_random)(a.engine_out, a.engine_log)

    # ── gates ───────────────────────────────────────────────────────────────
    errors = []
    if eng['engine_cases'] != coding['cases'] or eng['engine_controls'] != coding['controls']:
        errors.append(f"engine analysed {eng['engine_cases']} cases / {eng['engine_controls']} controls; "
                      f"phenotype file has {coding['cases']} / {coding['controls']}")
    if a.model == 'fixed' and eng['engine_covariates'] != a.n_covariates:
        errors.append(f"plink2 loaded {eng['engine_covariates']} covariates, expected {a.n_covariates}")
    if a.model == 'random':
        if eng['impute_method'] != 'mean':
            errors.append(f"SAIGE impute_method is {eng['impute_method']}, expected mean")
        if not eng['loco_chr6']:
            errors.append('SAIGE log does not confirm LOCO for chromosome 6')
    if not res['marker'].is_unique:
        errors.append('a marker has more than one result row')
    extra = set(res['marker']) - set(info['marker'])
    absent = set(info['marker']) - set(res['marker'])
    # An engine may skip a marker with no minor allele in the samples it analysed
    # (possible only in the case-only set, which has no controls). Any other
    # absent marker is an error.
    mono = set(info.loc[info['mac'] == 0, 'marker'])
    if extra:
        errors.append(f'{len(extra)} engine markers are not in marker_info')
    if absent - mono:
        errors.append(f'{len(absent - mono)} polymorphic markers have no engine result')
    if set(tests['marker']) != set(info['marker']):
        errors.append('the tests table and marker_info list different markers')
    if errors:
        sys.exit('ERROR: ' + '\n  '.join(errors))

    # ── one schema ──────────────────────────────────────────────────────────
    df = info.merge(res, on='marker', how='left', validate='one_to_one')
    df.loc[df['marker'].isin(absent), 'errcode'] = 'MONOMORPHIC'
    df.insert(3, 'model', a.model)
    df = df.merge(tests[['marker', 'test_id', 'mirror_of']], on='marker', how='left', validate='one_to_one')
    p = df['p'].to_numpy(dtype=float)
    for t in V_THRESHOLDS:
        df[f'sig_{t}'] = np.nan_to_num(p, nan=1.0) < float(th[f'p_{t}'])
    df = df.sort_values('bp').reset_index(drop=True)
    df.to_csv(a.out_assoc, sep='\t', index=False, float_format='%.6g')

    # ── scan summary ────────────────────────────────────────────────────────
    ok = ~np.isnan(p)
    errs = df.loc[df.errcode != '.', 'errcode'].value_counts().to_dict()
    top = df.loc[df['p'].idxmin()] if ok.any() else None
    qc = {'case_set': a.case_set, 'cohort': a.cohort, 'class': a.cls, 'model': a.model,
          'samples': coding['samples'], 'cases': coding['cases'], 'controls': coding['controls'],
          'markers_tested': m,
          'positions': int(df[['gene', 'position']].drop_duplicates().shape[0]) if a.cls == 'residue' else np.nan,
          'with_p': int(ok.sum()),
          'errcodes': ';'.join(f'{k}:{v}' for k, v in sorted(errs.items())) or '.',
          'lambda_gc': lambda_gc(p),
          'tests': int(th['tests']), 'mirror_pairs': int(th['mirror_pairs']),
          'me_gec': float(th['me_gec']),
          **{f'p_{t}': float(th[f'p_{t}']) for t in V_THRESHOLDS},
          **{f'n_{t}': int(df[f'sig_{t}'].sum()) for t in V_THRESHOLDS},
          **{f'tests_{t}': int(df.loc[df[f'sig_{t}'], 'test_id'].nunique()) for t in V_THRESHOLDS},
          'min_p': float(np.nanmin(p)) if ok.any() else np.nan,
          'top_marker': top['marker'] if top is not None else '',
          'top_or': float(top['or']) if top is not None else np.nan,
          'engine_cases': eng['engine_cases'], 'engine_controls': eng['engine_controls'],
          'engine_covariates': eng['engine_covariates'] if a.model == 'fixed' else a.n_covariates,
          'impute_method': eng.get('impute_method', 'none (complete case per marker)'),
          'variance_ratio': eng.get('variance_ratio', np.nan),
          'n_spa': eng.get('n_spa', np.nan)}
    pd.DataFrame([qc]).to_csv(a.out_qc, sep='\t', index=False, float_format='%.6g')
    print(f"[collect] {a.case_set} {a.cohort} {a.cls} {a.model}: {m} markers, {qc['tests']} tests, "
          f"{qc['with_p']} with P, lambda_GC {qc['lambda_gc']:.3f}, tests significant GEC / max-T / Bonferroni "
          f"{qc['tests_gec']} / {qc['tests_maxt']} / {qc['tests_bonf']}, min P {qc['min_p']:.3g} ({qc['top_marker']})")


if __name__ == '__main__':
    main()
