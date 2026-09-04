#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : ONE y AXIS PER VARIANT SET, decided once and written down.
#
#           The three cohorts of a variant set are the same statistic on nested
#           data, so a reader compares point heights across their three scan
#           figures. That only works if the three share a y axis, and a figure
#           drawn from its own cell cannot know the other two. This step reads
#           every cell of the scan and records, per variant set, the largest
#           -log10 P any cohort reached under either statistic.
#
#           IT SHARES THE DATA MAXIMUM, NOT THE AXIS. Each figure still applies
#           its own Bonferroni floor and draws its own threshold lines: the
#           denominators differ between cohorts, so the lines sit at different
#           heights and the figures are NOT interchangeable. What is shared is
#           the scale a point is measured against.
# Component: assoc_gene
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
import vocab as V                      # noqa: E402

NEEDED = ('cohort', 'stratum', 'method', 'pvalue')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gene-scan', nargs='+', required=True,
                   help='every SCAN gene_scan.tsv; the cohort and stratum are columns of the '
                        'file, so their staged names carry no meaning')
    p.add_argument('--out-tsv', default='scan_scale.tsv')
    return p.parse_args()


def main():
    a = parse_args()
    frames = []
    for f in a.gene_scan:
        d = pd.read_csv(f, sep='\t', usecols=lambda c: c in NEEDED,
                        dtype={'cohort': str, 'stratum': str, 'method': str})
        missing = [c for c in NEEDED if c not in d.columns]
        if missing:
            raise SystemExit(f'ABORT: {f} has no column(s) {missing}')
        frames.append(d)
    if not frames:
        raise SystemExit('ABORT: no gene_scan.tsv given')
    scan = pd.concat(frames, ignore_index=True)
    scan['nlp'] = V.nlp(pd.to_numeric(scan['pvalue'], errors='coerce'))
    scan = scan[np.isfinite(scan['nlp'])]
    if not len(scan):
        raise SystemExit('ABORT: no usable p-value in any gene_scan.tsv')

    rows = []
    for stratum, d in scan.groupby('stratum', sort=True):
        top = d.loc[d['nlp'].idxmax()]
        rows.append({'stratum': stratum,
                     'y_data_max': f'{float(top["nlp"]):.10g}',
                     'source_cohort': top['cohort'],
                     'source_method': top['method'],
                     'n_cohorts': d['cohort'].nunique(),
                     'n_methods': d['method'].nunique(),
                     'n_tests': len(d)})
    out = pd.DataFrame(rows)
    out.to_csv(a.out_tsv, sep='\t', index=False)
    for r in rows:
        print(f'[scan_scale] {r["stratum"]}: y_data_max {r["y_data_max"]} '
              f'({r["source_cohort"]}, {r["source_method"]}) over {r["n_cohorts"]} cohort(s) '
              f'x {r["n_methods"]} statistic(s), {r["n_tests"]:,} tests')
    print(f'[scan_scale] {len(rows)} variant set(s) -> {a.out_tsv}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:                       # noqa: BLE001
        print(f'Error in scan_scale: {e}', file=sys.stderr)
        sys.exit(1)
