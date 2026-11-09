#!/usr/bin/env python3
"""Platform sensitivity, the summary tables across cohorts (reads results/; run in the report stage).

  platform.signals.tsv    GWAS leads and DQB1:57:V: pooled, per-stratum OR, Q / P / empirical P / I2, FE and RE
                          estimates, LOPO (z observed vs expected, difference test), case-only P, per cohort
  platform.genes.tsv      Tier 1-2 genes x variant set x method: P per subset, carriers by platform, Firth carrier
                          OR (descriptive), carriers x platform permutation P
  platform.positions.tsv  residue positions: omnibus P per subset (pooled, strata, LOPO) and per cohort
  platform.null_fits.tsv  per cohort x subset: cases, controls, tau, variance ratio of the SAIGE null model
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402
import report_lib as RL  # noqa: E402


def last_tau(log):
    txt = Path(log).read_text(errors='replace')
    m = re.findall(r'Tau:\s*\n\[1\]\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)', txt)
    return float(m[-1][1]) if m else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results', required=True)
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    res = Path(a.results)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    leads = RL.gwas_leads(RL.read_signals(res), alt=True)
    targets = [(r.unit, r.role, r.variant_id) for r in leads.itertuples()] + [(RL.V57, 'HLA residue', RL.V57)]

    rows, fits = [], []
    for c in RL.cohorts_with(res, lambda r, c: RL.platform_table(r, c, 'heterogeneity')):
        het = pd.read_csv(RL.platform_table(res, c, 'heterogeneity'), sep='\t', dtype={'variant_id': str}).set_index('variant_id')
        st = pd.read_csv(RL.platform_table(res, c, 'strata'), sep='\t', dtype={'variant_id': str})
        lo = pd.read_csv(RL.platform_table(res, c, 'lopo'), sep='\t', dtype={'variant_id': str})
        co = pd.read_csv(RL.platform_table(res, c, 'case_only'), sep='\t', dtype={'marker': str}).set_index('marker')
        for unit, role, vid in targets:
            if vid not in het.index:
                continue
            h = het.loc[vid]
            row = {'cohort': c, 'signal': unit, 'role': role, 'marker': vid,
                   'or_pooled': np.exp(h.beta_pooled), 'p_pooled': h.p_pooled,
                   'Q': h.get('Q'), 'p_Q': h.get('p_Q'), 'p_Q_empirical': h.get('p_Q_empirical'),
                   'I2': h.get('I2'), 'I2_l95': h.get('I2_l95'), 'I2_u95': h.get('I2_u95'),
                   'or_fe': np.exp(h.get('beta_fe')), 'or_re': np.exp(h.get('beta_re')),
                   'or_re_l95': np.exp(h.get('re_l95')), 'or_re_u95': np.exp(h.get('re_u95')), 'tau2': h.get('tau2'),
                   'p_case_only': co.loc[vid, 'p_lrt'] if vid in co.index else np.nan}
            for sname in RL.PLATFORM_ONLY:
                x = st[(st.variant_id == vid) & (st.subset == sname)]
                if len(x):
                    x = x.iloc[0]
                    row[f'or_{sname}'] = x['or']
                    row[f'or_{sname}_ci'] = f"{x.or_l95:.2f}-{x.or_u95:.2f}"
                    row[f'p_{sname}'] = x.p
            for _, l in lo[lo.variant_id == vid].iterrows():
                k = l.left_out
                row[f'or_without_{k}'] = np.exp(l.beta)
                row[f'z_without_{k}'] = l.z_obs
                row[f'z_expected_without_{k}'] = l.z_expected
                row[f'p_without_{k}'] = l.p
                row[f'p_diff_{k}'] = l.get('p_diff', np.nan)
            rows.append(row)
        for f in sorted((res / '08.subsets' / c / 'null').glob('*/null.varianceRatio.txt')):
            s = f.parent.name
            cnt = pd.read_csv(RL.subset_counts_file(res, c, s), sep='\t').set_index('subset')
            fits.append({'cohort': c, 'subset': s, 'cases': int(cnt.loc[s, 'cases_random']),
                         'controls': int(cnt.loc[s, 'controls_random']),
                         'variance_ratio': float(f.read_text().split()[0]), 'tau': last_tau(f.parent / 'null.fit.log')})
    pd.DataFrame(rows).to_csv(out / 'platform.signals.tsv', sep='\t', index=False)
    pd.DataFrame(fits).to_csv(out / 'platform.null_fits.tsv', sep='\t', index=False)

    g = []
    for c in K.COHORTS:
        for f in sorted((res / '09.platform_assoc' / c).glob('gene.*.tsv')):
            if f.name.startswith('gene.ctrl'):          # control sources: control_report.py
                continue
            x = pd.read_csv(f, sep='\t')
            x['cohort'] = c
            g.append(x)
    if g:
        g = pd.concat(g)
        keep = ['cohort', 'subset', 'gene', 'stratum', 'method', 'n_var', 'n_informative', 'pvalue', 'or', 'or_l95',
                'or_u95', 'carriers_case', 'carriers_control']
        gt = g[[k for k in keep if k in g.columns]]
        fir = [pd.read_csv(res / '09.platform_assoc' / c / 'carriers_firth.tsv', sep='\t')
               for c in RL.cohorts_with(res, lambda r, c: r / '09.platform_assoc' / c / 'carriers_firth.tsv')]
        if fir:
            fir = pd.concat(fir)[['cohort', 'subset', 'gene', 'stratum', 'or_firth', 'l95', 'u95', 'p_firth']]
            gt = gt.merge(fir.rename(columns={'l95': 'or_firth_l95', 'u95': 'or_firth_u95'}),
                          on=['cohort', 'subset', 'gene', 'stratum'], how='left')
        perm = []
        for c in RL.cohorts_with(res, lambda r, c: RL.platform_table(r, c, 'case_only')):
            x = pd.read_csv(RL.platform_table(res, c, 'case_only'), sep='\t', dtype={'marker': str})
            x = x[x.marker.str.contains('__')]
            x[['gene', 'stratum']] = x.marker.str.split('__', expand=True)
            x['cohort'] = c
            perm.append(x[['cohort', 'gene', 'stratum', 'carriers_by_platform', 'p_perm', 'p_lrt']])
        if perm:
            gt = gt.merge(pd.concat(perm), on=['cohort', 'gene', 'stratum'], how='left')
        gt.to_csv(out / 'platform.genes.tsv', sep='\t', index=False)

    om, w = RL.omnibus_wide([(f.name.split('.')[1], f)
                             for f in sorted((res / '09.platform_assoc' / '_omnibus').glob('omnibus.*.tsv'))])
    if om is not None:
        pooled = om[om.subset == 'pooled'].set_index(['cohort', 'unit'])
        w['threshold'] = [float(pooled.loc[(r.cohort, r.unit), 'threshold']) if (r.cohort, r.unit) in pooled.index
                          else np.nan for r in w.itertuples()]
        w.to_csv(out / 'platform.positions.tsv', sep='\t', index=False)
    print(f'[platform_report] {len(rows)} signal rows, {len(fits)} null fits -> {out}')


if __name__ == '__main__':
    main()
