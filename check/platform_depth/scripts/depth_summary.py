#!/usr/bin/env python3
"""Stage 01 — depth per genome and platform, the down-sampling target and fractions, and the
read-length audit (regenerated: the info sheet's Observed_Depth for DNBSeq-T7 had assumed
150 bp reads; Edit 2 of the sheet corrected it).

Target D* = the median effective depth (probe windows, duplicates excluded) of the full-mainland
HiSeq X controls. Each deep case (NovaSeq, DNBSeq-G400RS 30x) keeps f_i = D* / D_i of its read
pairs (user decision 2026-10-30); every other genome is untouched.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402


def region_len(reg):
    _, span = reg.split(':')
    s, e = (int(x) for x in span.split('-'))
    return e - s + 1


def wmean(g):
    w = g.region.map(region_len)
    return float(np.average(g.meandepth.astype(float), weights=w))


def q(x, p):
    return float(np.percentile(x, p))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--depth', nargs='+', required=True, help='depth_probe.py chunk outputs')
    ap.add_argument('--samples', required=True, help='sample_table.tsv')
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    d = pd.concat([K.read_tsv(f) for f in a.depth], ignore_index=True)
    st = K.read_tsv(a.samples)
    missing = set(st.sample_id) - set(d.sample_id)
    if missing:
        sys.exit(f'ABORT: no depth for {len(missing)} genomes, e.g. {sorted(missing)[:5]}')
    n_reg = d.groupby(['sample_id', 'set', 'duplicates']).size().unstack(['set', 'duplicates'])
    if n_reg.nunique().max() != 1:
        sys.exit('ABORT: genomes differ in the number of regions measured')

    per = []
    for sid, g in d.groupby('sample_id'):
        pr = g[g.set == 'probe']
        per.append({'sample_id': sid,
                    'probe_depth': wmean(pr[pr.duplicates == 'excluded']),
                    'probe_depth_with_dup': wmean(pr[pr.duplicates == 'kept']),
                    'target_depth': wmean(g[g.set == 'target']),
                    'readlen_median': float(g.readlen_median.iloc[0]), 'readlen_max': int(g.readlen_max.iloc[0])})
    per = pd.DataFrame(per)
    per['dup_share'] = 1 - per.probe_depth / per.probe_depth_with_dup
    per = st[['sample_id', 'status', 'platform', 'deep', 'observed_depth', 'batch'] +
             [f'in_{c}' for c in K.COHORTS]].merge(per, on='sample_id')
    per['observed_depth'] = pd.to_numeric(per.observed_depth, errors='coerce')
    # the audit: the read length the sheet's depth implies, L_real x Observed / CRAM (duplicates kept)
    per['readlen_implied'] = per.readlen_median * per.observed_depth / per.probe_depth_with_dup
    K.write_tsv(per.round(4), out / 'depth_per_sample.tsv')

    full = per[per.in_full_mainland == '1']
    ctrl = full[full.platform == K.CONTROL_PLATFORM]
    target = float(np.median(ctrl.probe_depth))

    rows = []
    for c in K.COHORTS:
        dc = per[per[f'in_{c}'] == '1']
        for pf in [K.CONTROL_PLATFORM] + K.CASE_PLATFORMS:
            g = dc[dc.platform == pf]
            rows.append({'cohort': c, 'platform': pf, 'status': 'control' if pf == K.CONTROL_PLATFORM else 'case',
                         'n': len(g), 'readlen_median': float(np.median(g.readlen_median)),
                         'depth_median': float(np.median(g.probe_depth)),
                         'depth_q25': q(g.probe_depth, 25), 'depth_q75': q(g.probe_depth, 75),
                         'depth_with_dup_median': float(np.median(g.probe_depth_with_dup)),
                         'dup_share_median': float(np.median(g.dup_share)),
                         'observed_depth_median': float(np.median(g.observed_depth)),
                         'observed_over_measured_median': float(np.median(g.observed_depth / g.probe_depth_with_dup)),
                         'readlen_implied_median': float(np.median(g.readlen_implied)),
                         'target_depth_median': float(np.median(g.target_depth))})
    plat = pd.DataFrame(rows)
    K.write_tsv(plat.round(4), out / 'platform_depth.tsv')

    deep = full[full.deep == '1'].copy()
    deep['fraction'] = target / deep.probe_depth
    if not ((deep.fraction > 0) & (deep.fraction < 1)).all():
        bad = deep[(deep.fraction <= 0) | (deep.fraction >= 1)]
        sys.exit(f'ABORT: fraction outside (0, 1) for {len(bad)} deep cases: {bad.sample_id.tolist()[:5]}')
    deep['expected_depth'] = deep.probe_depth * deep.fraction
    fr = deep[['sample_id', 'platform', 'batch', 'probe_depth', 'probe_depth_with_dup', 'observed_depth',
               'fraction', 'expected_depth']]
    K.write_tsv(fr.round(6), out / 'fractions.tsv')

    # local depth at every calling interval: platform medians, and the deep cases after scaling
    loc = d[d.set == 'target'].merge(per[['sample_id', 'platform', 'in_full_mainland']], on='sample_id')
    loc = loc[loc.in_full_mainland == '1'].copy()
    loc['meandepth'] = loc.meandepth.astype(float)
    loc = loc.merge(fr[['sample_id', 'fraction']], on='sample_id', how='left')
    loc['after'] = loc.meandepth * loc.fraction.fillna(1.0)
    lt = loc.pivot_table(index='region', columns='platform', values='meandepth', aggfunc='median')
    lt.columns = [f'median_{K.PLATFORM_SLUG[c]}' for c in lt.columns]
    da = loc[loc.fraction.notna()].groupby('region').after.median().rename('median_deep_after')
    lt = lt.join(da)
    lt['deep_after_over_controls'] = lt.median_deep_after / lt.median_hiseqx
    K.write_tsv(lt.reset_index().round(4), out / 'local_depth.tsv')

    summ = {'target_depth': round(target, 4), 'target_source': 'median probe depth, duplicates excluded, '
            f'full-mainland {K.CONTROL_PLATFORM} controls (n = {len(ctrl)})',
            'deep_cases': len(deep), 'fraction_median': round(float(deep.fraction.median()), 4),
            'fraction_by_platform': {pf: [round(float(g.fraction.min()), 4), round(float(g.fraction.median()), 4),
                                          round(float(g.fraction.max()), 4)] for pf, g in deep.groupby('platform')},
            'local_ratio_deep_after_over_controls': [round(float(lt.deep_after_over_controls.min()), 3),
                                                     round(float(lt.deep_after_over_controls.median()), 3),
                                                     round(float(lt.deep_after_over_controls.max()), 3)]}
    (out / 'depth_summary.json').write_text(json.dumps(summ, indent=2) + '\n')
    print(json.dumps(summ, indent=2))
    print(plat[plat.cohort == 'full_mainland'].round(3).to_string(index=False))


if __name__ == '__main__':
    main()
