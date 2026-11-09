#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The DVT status of every sample, from the v6 sample sheet, in the one
#           three-level coding the case sets use: yes / no / unknown.
#
#           DVT_Presence is recorded for cases only ('Yes', 'No', 'Unknown' or
#           blank). 'Unknown' and blank both become 'unknown'; any other value
#           stops the run rather than being guessed at. Controls carry no DVT
#           record and are never split: every case set uses all of them.
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import sys

import pandas as pd

MAP = {'yes': 'yes', 'no': 'no', 'unknown': 'unknown', '': 'unknown'}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--info', required=True, help='cteph_agp3k.v6 sample sheet (.xlsx)')
    p.add_argument('--id-col', default='ID_JHRPv6')
    p.add_argument('--dvt-col', default='DVT_Presence')
    p.add_argument('--out', required=True)
    return p.parse_args()


def main():
    a = parse_args()
    df = pd.read_excel(a.info, dtype=str)
    for c in (a.id_col, a.dvt_col):
        if c not in df.columns:
            sys.exit(f'ERROR: column {c} not in {a.info}')
    df = df[df[a.id_col].notna()]
    if not df[a.id_col].is_unique:
        sys.exit(f'ERROR: {a.id_col} is not unique in {a.info}')
    raw = df[a.dvt_col].fillna('').str.strip()
    bad = sorted(set(raw.str.lower()) - set(MAP))
    if bad:
        sys.exit(f'ERROR: unexpected {a.dvt_col} values: {bad}')
    out = pd.DataFrame({'IID': df[a.id_col], 'dvt': raw.str.lower().map(MAP), 'dvt_raw': raw})
    out.to_csv(a.out, sep='\t', index=False)
    print(f'[case_sets] {len(out)} samples; DVT: ' +
          ', '.join(f'{k} {v}' for k, v in out.dvt.value_counts().items()))


if __name__ == '__main__':
    main()
