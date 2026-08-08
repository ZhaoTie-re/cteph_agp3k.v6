#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Recommend a minAC threshold from the merged per-minAC QC table.
#           The depth-(and PC-)adjusted apparent phenotype effect beta_group is
#           V-SHAPED in minAC: inflated at minAC 0-1, dipping to a NON-significant
#           trough (the decoupling point), then climbing again as higher-count
#           variants re-introduce platform inflation. The recommendation is the
#           trough:
#             = the SMALLEST minAC (>= floor) at which beta_group is non-significant
#               (p >= alpha AND 0 in the 95% CI).
#           This is the left edge of the non-significant band — it maximises power
#           (most variants) while decoupling the burden from depth/platform.
#           If nothing is non-significant, fall back to argmin|beta_group| (the
#           least-confounded threshold) and flag calibrated=false.
#
#           PRIMARY model = no-PC. beta_group here is a DETECTOR for the technical
#           artifact, and with zero platform overlap that artifact IS the case/control
#           contrast — so a covariate carrying that contrast hides part of what the
#           detector is for. Measured, PC adjustment moves the reading by 1.9-2.6x its
#           own size and in cohort-dependent directions, so PC does not drive the
#           threshold; the +PC model is reported as a SENSITIVITY check and we flag
#           whether the two agree at the recommended minAC. See docs/METHODS.md S5.4
#           and results/_comparison/decision_axis.png.
# Project : cteph_agp3k.v6 / tuning.rv  (rare-variant depth-confounding QC)
# Used by : tuning.rv.nf  process SELECT_MINAC
# ---------------------------------------------------------------------------
import argparse
import sys

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description='Recommend the V-trough (decoupling) minAC.')
    p.add_argument('--qc-stats', required=True, help='Merged per-minAC QC table (all_qc_summary_stats.tsv).')
    p.add_argument('--out-tsv', required=True, help='Recommendation sidecar (key<TAB>value).')
    p.add_argument('--alpha', type=float, default=0.05, help='Non-significance level (default 0.05).')
    p.add_argument('--min-floor', type=int, default=1,
                   help='Smallest minAC eligible to be recommended (default 1; excludes singletons at 0).')
    p.add_argument('--prefer', choices=['nopc', 'pc'], default='nopc',
                   help='PRIMARY model for the recommendation (default nopc; pc reported as sensitivity).')
    return p.parse_args()


def _cols(model):
    s = 'pc' if model == 'pc' else 'noPC'
    return f'Beta_group_{s}', f'Beta_group_{s}_L', f'Beta_group_{s}_U', f'P_group_{s}'


def _analyze(elig, model, alpha):
    """Return the trough recommendation for one model, or None if that model is absent."""
    b, l, u, p = _cols(model)
    if b not in elig.columns or not elig[b].notna().any():
        return None
    d = elig.dropna(subset=[b, l, u, p])
    if d.empty:
        return None
    nonsig = (d[p] >= alpha) & (d[l] <= 0) & (d[u] >= 0)      # 0 inside the 95% CI, non-significant
    band = d.loc[nonsig, 'MinAC_Threshold'].astype(int).tolist()
    trough = int(d.loc[d[b].abs().idxmin(), 'MinAC_Threshold'])  # least-confounded threshold
    if band:
        rec, calibrated = min(band), True                     # left edge of the non-sig trough = max power
    else:
        rec, calibrated = trough, False                       # fallback: least-confounded (flag)
    row = d[d['MinAC_Threshold'] == rec].iloc[0]
    return {
        'model': 'pc' if model == 'pc' else 'noPC', 'rec': rec, 'calibrated': calibrated,
        'band': band, 'trough': trough,
        'beta': float(row[b]), 'lo': float(row[l]), 'hi': float(row[u]), 'p': float(row[p]),
        'nvar': int(row['Variant_Count']) if ('Variant_Count' in row and pd.notna(row['Variant_Count'])) else -1,
    }


def _nonsig_at(elig, model, minac, alpha):
    """Is `minac` non-significant (0 in CI) under `model`? None if unknown."""
    b, l, u, p = _cols(model)
    r = elig[elig['MinAC_Threshold'] == minac]
    if not len(r) or b not in elig.columns or pd.isna(r.iloc[0][p]):
        return None
    r = r.iloc[0]
    return bool((r[p] >= alpha) and (r[l] <= 0) and (r[u] >= 0))


def main():
    args = parse_args()
    df = pd.read_csv(args.qc_stats, sep='\t').sort_values('MinAC_Threshold').reset_index(drop=True)
    if df.empty:
        raise ValueError('Empty qc-stats table.')
    elig = df[df['MinAC_Threshold'] >= args.min_floor].reset_index(drop=True)
    pc_source = str(df['PC_source'].dropna().iloc[0]) if ('PC_source' in df.columns
                                                          and df['PC_source'].notna().any()) else 'none'

    prefer = args.prefer
    res = _analyze(elig, prefer, args.alpha)
    if res is None and prefer == 'pc':                 # PCs unavailable -> fall back to no-PC
        prefer = 'nopc'
        res = _analyze(elig, prefer, args.alpha)
    if res is None:
        raise ValueError('No usable beta_group columns in the qc-stats table.')

    other = 'pc' if prefer == 'nopc' else 'nopc'
    res_o = _analyze(elig, other, args.alpha)
    agree = _nonsig_at(elig, other, res['rec'], args.alpha) if res_o is not None else None

    out = {
        'recommended_minac': res['rec'],
        'preferred_model': res['model'],
        'calibrated': str(res['calibrated']).lower(),
        'rule': f'smallest minAC>={args.min_floor} with p_group>={args.alpha:g} AND 0 in 95%CI (V-trough left edge)',
        'alpha': f'{args.alpha:g}', 'min_floor': args.min_floor, 'pc_source': pc_source,
        'trough_minac': res['trough'],
        'nonsig_band': ','.join(map(str, res['band'])) if res['band'] else 'none',
        'beta_group_at_rec': f'{res["beta"]:.6g}',
        'ci_low_at_rec': f'{res["lo"]:.6g}', 'ci_high_at_rec': f'{res["hi"]:.6g}',
        'p_group_at_rec': f'{res["p"]:.6g}',
        'variant_count_at_rec': res['nvar'],
        # sensitivity (other model)
        'sensitivity_model': res_o['model'] if res_o else 'none',
        'sensitivity_rec': res_o['rec'] if res_o else -1,
        'sensitivity_nonsig_band': (','.join(map(str, res_o['band'])) if (res_o and res_o['band'])
                                    else ('none' if res_o else 'na')),
        'models_agree_at_rec': ('na' if agree is None else str(agree).lower()),
    }
    with open(args.out_tsv, 'w') as f:
        f.write('key\tvalue\n')
        for k, v in out.items():
            f.write(f'{k}\t{v}\n')

    tag = 'CALIBRATED (V-trough)' if res['calibrated'] else 'NOT-CALIBRATED (fell back to least-confounded minAC)'
    print(f'Recommended minAC : {res["rec"]}   [{tag}]   primary model = {res["model"]}')
    print(f'  non-sig band    : {out["nonsig_band"]}   (trough |beta| at minAC={res["trough"]})')
    print(f'  beta_group@rec  : {out["beta_group_at_rec"]}  '
          f'(95% CI {out["ci_low_at_rec"]}..{out["ci_high_at_rec"]}), p={out["p_group_at_rec"]}')
    if out['variant_count_at_rec'] >= 0:
        print(f'  variants kept   : {out["variant_count_at_rec"]:,}')
    if res_o is not None:
        agree_txt = {'true': 'AGREE', 'false': 'DISAGREE', 'na': 'n/a'}[out['models_agree_at_rec']]
        print(f'  sensitivity     : {other} model recommends minAC={res_o["rec"]}; '
              f'at minAC={res["rec"]} the two models {agree_txt}')
    print(f'  sidecar         : {args.out_tsv}')
    print(res['rec'])          # bare number on the last line for Nextflow capture


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in select_minac: {e}', file=sys.stderr)
        sys.exit(1)
