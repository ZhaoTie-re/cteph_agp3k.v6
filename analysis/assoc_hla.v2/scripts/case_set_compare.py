#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Put the four case sets side by side: all cases, cases with DVT and
#           cases without DVT (each against the same controls), and the
#           case-only test of DVT yes against DVT no.
#
#   scan_qc_case_sets.tsv  the 48 scan summaries (4 case sets x 3 cohorts x
#                          2 classes x 2 models)
#   dvt_contrast.tsv       one row per (class, cohort, model, marker): OR, P and
#                          significance at each threshold (GEC primary, max-T,
#                          conventional Bonferroni) in each case set. The DVT-yes
#                          and DVT-no columns share their controls, so their
#                          estimates are correlated and their gap is not a
#                          test; the dvt_case_only columns ARE the test of a
#                          difference
#   dvt_hits.tsv           the dvt_contrast rows of every marker significant
#                          at any threshold in a DVT set or in the case-only
#                          test, in any cohort or model
#
#           Tables only: the figures are drawn by figures.py (02.signals,
#           05.dvt).
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import pandas as pd                      # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import vocab as V                        # noqa: E402

SIG = [f'sig_{t}' for t in V.THRESHOLDS]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assoc', nargs='+', required=True, help="each case set's assoc_all.tsv")
    p.add_argument('--qc', nargs='+', required=True, help="each case set's scan_qc_all.tsv")
    p.add_argument('--out-dir', default='.')
    return p.parse_args()


def main():
    a = parse_args()
    out = Path(a.out_dir)
    assoc = pd.concat([pd.read_csv(f, sep='\t', dtype={'marker': str}) for f in a.assoc], ignore_index=True)
    qc = pd.concat([pd.read_csv(f, sep='\t') for f in a.qc], ignore_index=True)
    if set(assoc.case_set) != set(V.CASE_SETS):
        sys.exit(f'ERROR: case sets present {sorted(set(assoc.case_set))}, expected {V.CASE_SETS}')
    co = {c: i for i, c in enumerate(V.COHORTS)}
    cs = {c: i for i, c in enumerate(V.CASE_SETS)}
    qc = qc.sort_values(['case_set', 'class', 'cohort', 'model'],
                        key=lambda s: s.map(cs) if s.name == 'case_set' else (s.map(co) if s.name == 'cohort' else s))
    qc.to_csv(out / 'scan_qc_case_sets.tsv', sep='\t', index=False, float_format='%.6g')

    key = ['class', 'cohort', 'model', 'marker']
    wide = None
    for c in V.CASE_SETS:
        d = assoc[assoc.case_set == c][key + ['gene', 'position', 'residue', 'n_case', 'freq_case',
                                              'freq_control', 'or', 'l95', 'u95', 'p'] + SIG]
        # Marker annotation and the control frequency come from 'all' once; the
        # case-only set's "control" frequency (DVT no) is its freq_control_dvt_case_only.
        shared = ['gene', 'position', 'residue', 'freq_control']
        if c != 'all':
            d = d.drop(columns=shared[:3] + ([] if c == 'dvt_case_only' else ['freq_control']))
        d = d.rename(columns={k: f'{k}_{c}' for k in d.columns if k not in key + (shared if c == 'all' else [])})
        wide = d if wide is None else wide.merge(d, on=key, how='outer')
    wide['position'] = wide['position'].astype('Int64')
    for c in V.CASE_SETS:
        wide[[f'{k}_{c}' for k in SIG]] = wide[[f'{k}_{c}' for k in SIG]].fillna(False).astype(bool)
    anysig = lambda c: wide[[f'{k}_{c}' for k in SIG]].any(axis=1)          # noqa: E731
    wide['called_dvt'] = anysig('dvt_yes') | anysig('dvt_no')
    wide['called_case_only'] = anysig('dvt_case_only')
    wide['min_p_dvt'] = wide[['p_dvt_yes', 'p_dvt_no']].min(axis=1)
    wide = wide.sort_values(['class', 'min_p_dvt']).reset_index(drop=True)
    wide.to_csv(out / 'dvt_contrast.tsv', sep='\t', index=False, float_format='%.6g')
    called = set(wide.loc[wide.called_dvt | wide.called_case_only, 'marker'])
    called_co = set(wide.loc[wide.called_case_only, 'marker'])
    hits = wide[wide.marker.isin(called)]
    hits.to_csv(out / 'dvt_hits.tsv', sep='\t', index=False, float_format='%.6g')

    print(f'[case_set_compare] {len(qc)} scans; {len(called)} markers significant in a DVT set or case-only '
          f'({len(called_co)} in the case-only test)')


if __name__ == '__main__':
    main()
