#!/usr/bin/env python3
"""Gates of the down-sampling arms, from their own outputs.

  G6  re-processing fidelity: concordance >= 99.5 % for baseline vs original, original GVCFs re-genotyped
      vs original, and baseline vs original GVCFs (deep cases and fidelity controls), and the call rate of
      the baseline within 0.5 percentage points of the original's
  G8  achieved depth: over the calling intervals, depth after / (f x depth before) has a median in
      0.95-1.05 and every genome in 0.85-1.15 (the subsampling did what was asked)
  G9  HLA re-typing at f = 1 reproduces the archived typing (every call identical), or the differences are
      listed and the f = 1 re-typing serves as the baseline
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--concordance', required=True, help='concordance.summary.tsv')
    ap.add_argument('--depth', nargs='+', required=True, help='SUBSAMPLE *.depth.tsv')
    ap.add_argument('--hla', required=True, help='hla.calls_long.tsv')
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    out = Path(a.out_dir)

    c = pd.read_csv(a.concordance, sep='\t')
    c = c[(c.kind == 'all') & (c.maf_class == 'all') & c.platform.isna()]
    fid = c[c.comparison.isin(['recall_vs_original', 'origgvcf_vs_original', 'recall_vs_origgvcf'])]
    rows = fid[['comparison', 'role', 'arm', 'samples', 'both_called', 'concordance', 'nonref_concordance',
                'callrate_ref', 'callrate_test']].to_dict('records')
    ok6 = bool((fid.concordance >= 0.995).all() and ((fid.callrate_test - fid.callrate_ref).abs() <= 0.005).all())
    json.dump({'gate': 'G6', 'rows': rows, 'pass': ok6}, open(out / 'G6.json', 'w'), indent=2, default=float)

    d = pd.concat([pd.read_csv(f, sep='\t') for f in a.depth], ignore_index=True)
    d = d[d.arm != 'baseline'].copy()
    d['ratio'] = d.depth_after / (d.fraction * d.depth_before)
    per = d.groupby(['sample_id', 'arm']).ratio.median()
    ok8 = bool(0.95 <= per.median() <= 1.05 and per.between(0.85, 1.15).all())
    json.dump({'gate': 'G8', 'genome_arms': int(len(per)), 'ratio_median': float(per.median()),
               'ratio_min': float(per.min()), 'ratio_max': float(per.max()), 'pass': ok8},
              open(out / 'G8.json', 'w'), indent=2)

    h = pd.read_csv(a.hla, sep='\t', dtype=str)
    b = h[h.comparison == 'baseline_vs_archive']
    diff = b[b.identical != '1']
    ok9 = bool(len(b) > 0 and len(diff) == 0)
    json.dump({'gate': 'G9', 'calls': int(len(b)), 'differing': int(len(diff)),
               'differences': diff[['sample_id', 'gene', 'state_ref', 'state_test', 'alleles_ref',
                                    'alleles_test']].to_dict('records'),
               'pass': ok9, 'note': 'if differing > 0, the f = 1 re-typing is the baseline for the depth effect'},
              open(out / 'G9.json', 'w'), indent=2)
    print(f'G6 {ok6}  G8 {ok8} (median {per.median():.4f}, range {per.min():.4f}-{per.max():.4f})  '
          f'G9 {ok9} ({len(diff)} of {len(b)} calls differ)')


if __name__ == '__main__':
    main()
