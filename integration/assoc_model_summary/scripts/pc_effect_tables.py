#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The evidence on what the PCs do on top of the full GRM, as tables
#           (the argument itself is to be agreed; this only states facts):
#             peak_changes.tsv       per cohort, the peaks called under only one
#                                    of the two models (no peak of the other
#                                    within the window), with gene and size
#             lead_or_stability.tsv  OR and P of the cross-cohort leads under
#                                    both models
#             pc_effect.json         per cohort: lambda_GC, variants at
#                                    P < 1e-5, peaks only in one model, the PC /
#                                    case-status test, OR shift summary
# Inputs  : stats.json of both models; assoc_saige model_peaks_annotation.tsv
#           per cohort and _comparison/tables/lead_bymodel.tsv;
#           pc_case_association.tsv
# Outputs : peak_changes.tsv, lead_or_stability.tsv, pc_effect.json
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohorts', required=True)
    ap.add_argument('--with-pc', default='grm_pc')
    ap.add_argument('--without-pc', default='grm_only')
    ap.add_argument('--peaks', required=True, help='assoc_saige model_peaks_annotation.tsv per cohort, comma-separated')
    ap.add_argument('--upstream-models', required=True)
    ap.add_argument('--lead-bymodel', required=True)
    ap.add_argument('--pc-assoc', required=True)
    ap.add_argument('--window-bp', type=int, default=500000)
    ap.add_argument('--out-peaks', required=True)
    ap.add_argument('--out-or', required=True)
    ap.add_argument('--out-json', required=True)
    a = ap.parse_args()
    cohorts = C.csv_list(a.cohorts)
    up = C.kv(a.upstream_models)
    mp, mn = a.with_pc, a.without_pc
    pc = pd.read_csv(a.pc_assoc, sep='\t').set_index('cohort')
    rows, summary = [], {}
    for c, f in zip(cohorts, C.csv_list(a.peaks)):
        d = pd.read_csv(f, sep='\t', dtype={'chrom': str})
        d['pos'] = pd.to_numeric(d.pos)
        both = {m: d[d.model == up[m]] for m in (mp, mn)}
        only = {}
        for m, other in ((mp, mn), (mn, mp)):
            o = both[other]
            keep = []
            for _, r in both[m].iterrows():
                near = o[(o.chrom == r.chrom) & ((o.pos - r.pos).abs() <= a.window_bp)]
                if not len(near):
                    keep.append(r)
                    rows.append(dict(cohort=c, only_in=m, only_in_label=C.MODELS[m]['short'], gene=r.Gene,
                                     chrom=r.chrom, pos=int(r.pos), variant_id=r.variant_id, rsID=r.rsID,
                                     tier=r.tier, n_sig_variants=int(r.n_sig_variants), P=float(r.P),
                                     OR=float(r.OR)))
            only[m] = keep
        sp, sn = (C.load_json(f'{c}.{m}.stats.json') for m in (mp, mn))
        summary[c] = dict(
            lambda_with_pc=sp['lambda_gc'], lambda_without_pc=sn['lambda_gc'],
            p_lt_1e5_with_pc=sp['n_p_lt_suggestive'], p_lt_1e5_without_pc=sn['n_p_lt_suggestive'],
            expected_p_lt_1e5=sp['expected_p_lt_suggestive'],
            peaks_with_pc=len(both[mp]), peaks_without_pc=len(both[mn]),
            peaks_only_with_pc=len(only[mp]), peaks_only_without_pc=len(only[mn]),
            sig_variants_in_peaks_only_without_pc=int(sum(r.n_sig_variants for r in only[mn])),
            largest_peak_only_without_pc=(max(only[mn], key=lambda r: r.n_sig_variants).Gene if only[mn] else ''),
            largest_peak_only_without_pc_n=(int(max(r.n_sig_variants for r in only[mn])) if only[mn] else 0),
            pc_case_lr_chi2=float(pc.loc[c, 'lr_chi2']), pc_case_lr_df=int(pc.loc[c, 'lr_df']),
            pc_case_lr_p=float(pc.loc[c, 'lr_p']), pc_case_n_pcs_p_lt_0_05=int(pc.loc[c, 'n_pcs_p_lt_0_05']))
    pd.DataFrame(rows).to_csv(a.out_peaks, sep='\t', index=False, float_format='%.6g')

    lb = pd.read_csv(a.lead_bymodel, sep='\t')
    lb = lb[lb.cohort.isin(cohorts)].copy()
    lb['dlnOR'] = np.log(lb.OR_main) - np.log(lb.OR_detect)
    lb['dnlp'] = -np.log10(lb.P_main) + np.log10(lb.P_detect)
    lb.to_csv(a.out_or, sep='\t', index=False, float_format='%.6g')
    for c in cohorts:
        s = lb[lb.cohort == c]
        summary[c].update(n_leads_compared=len(s), median_abs_dlnOR=float(s.dlnOR.abs().median()),
                          max_abs_dlnOR=float(s.dlnOR.abs().max()))
    with open(a.out_json, 'w') as fh:
        json.dump(dict(with_pc=mp, without_pc=mn, window_bp=a.window_bp, cohorts=summary), fh, indent=2)
    C.log(f'pc effect: {len(rows)} model-specific peaks; {len(lb)} leads compared')


if __name__ == '__main__':
    main()
