#!/usr/bin/env python3
"""The down-sampling fractions of every down-sampled case: the deep cases' (DEPTH_SUMMARY, unchanged) and the
DNBSeq-T7 cases' (user, 2026-11-03).

Every case platform but DNBseq-G400RS 15x is brought to the controls' depth, so that each platform's median
effective depth equals the HiSeq X controls' median D*:
  NovaSeq 30x, DNBSeq-G400RS 30x   each genome to D* (f_i = D* / D_i; fractions.tsv of DEPTH_SUMMARY, copied
                                   line for line, so the arms already run keep their fraction and their cache)
  DNBSeq-T7 30x                    one fraction for the platform, f = D* / (the full-mainland DNBSeq-T7 median),
                                   which keeps the platform's own spread of depth and puts its median at D*
  DNBseq-G400RS 15x                not touched (below D*; depth cannot be added)
The DNBSeq-T7 cases are those of full mainland, the cohort that holds the other two.

Output: one table in the layout of fractions.tsv (the deep cases first, then the DNBSeq-T7 cases by sample ID).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

T7 = 'DNBSeq-T7 30x'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fractions', required=True, help='DEPTH_SUMMARY fractions.tsv (the deep cases)')
    ap.add_argument('--per-sample', required=True, help='DEPTH_SUMMARY depth_per_sample.tsv')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    lines = Path(a.fractions).read_text().splitlines()
    fr = pd.read_csv(a.fractions, sep='\t', dtype=str)
    d_star = float(fr.expected_depth.astype(float).median())
    ps = pd.read_csv(a.per_sample, sep='\t', dtype={'sample_id': str})
    full = ps[ps.in_full_mainland == 1]
    ctrl = full[full.status == 'control'].probe_depth.median()
    if not np.isclose(ctrl, d_star, rtol=1e-5):
        sys.exit(f'ABORT: the controls\' median {ctrl} differs from D* {d_star} of fractions.tsv')
    t7 = full[(full.status == 'case') & (full.platform == T7)].sort_values('sample_id')
    if fr.sample_id.isin(t7.sample_id).any():
        sys.exit('ABORT: a DNBSeq-T7 case is already in fractions.tsv')
    f = round(d_star / t7.probe_depth.median(), 5)
    if not 0 < f < 1:
        sys.exit(f'ABORT: DNBSeq-T7 fraction {f} outside (0, 1)')
    cols = lines[0].split('\t')
    rows = []
    for r in t7.itertuples(index=False):
        v = {'sample_id': r.sample_id, 'platform': r.platform, 'batch': r.batch, 'probe_depth': f'{r.probe_depth:g}',
             'probe_depth_with_dup': f'{r.probe_depth_with_dup:g}', 'observed_depth': f'{r.observed_depth:g}',
             'fraction': f'{f:g}', 'expected_depth': f'{f * r.probe_depth:.6f}'}
        rows.append('\t'.join(v[c] for c in cols))
    Path(a.out).write_text('\n'.join(lines + rows) + '\n')
    print(f'[t7_fractions] D* {d_star:.4f}; DNBSeq-T7 median {t7.probe_depth.median():.4f} -> f = {f}; '
          f'{len(fr)} deep cases + {len(t7)} DNBSeq-T7 cases')


if __name__ == '__main__':
    main()
