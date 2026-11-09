#!/usr/bin/env python3
"""The down-sampling design: which genome is processed in which arm, with which fraction and seed.

Arms
  baseline   f = 1: every read, through the same new processing (the reference for the depth effect)
  r01..rNN   f_i = D* / D_i, seed = replicate number (the same seed for variant and HLA reads)
Roles
  deep       the 89 NovaSeq + DNBSeq-G400RS 30x cases: baseline + every replicate, swapped into the
             association filesets and re-typed for HLA
  fidelity   HiSeq X controls with a surviving original GVCF, baseline only: re-processing fidelity at
             control depth (never swapped)
Original GVCFs (orig_gvcfs.tsv) of deep cases and fidelity controls are joint-genotyped beside the
arms, so the joint-genotyping step itself is checked against the cohort's calls.

--mode pilot keeps 2 v5 NovaSeq cases with original GVCFs, 1 NovaSeq without, 1 DNBSeq-G400RS 30x
case and 2 fidelity controls (first by sample ID in each class: no choice by outcome).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--fractions', required=True)
    ap.add_argument('--samples', required=True)
    ap.add_argument('--mode', choices=['pilot', 'full'], required=True)
    ap.add_argument('--replicates', type=int, required=True)
    ap.add_argument('--fidelity-controls', type=int, default=20)
    ap.add_argument('--seed', type=int, default=20261030)
    ap.add_argument('--out-plan', required=True)
    ap.add_argument('--out-gvcfs', required=True)
    ap.add_argument('--out-deep', required=True)
    a = ap.parse_args()

    fr = K.read_tsv(a.fractions)
    st = K.read_tsv(a.samples).set_index('sample_id')
    deep = fr.merge(st[['cram', 'crai', 'orig_gvcf_dir', 'batch']].reset_index(), on='sample_id', suffixes=('', '_st'))
    ctrl = st[(st.status == 'control') & (st.in_full_mainland == '1') & (st.orig_gvcf_dir != '')].reset_index()
    if a.mode == 'pilot':
        nov = deep[deep.platform == 'NovaSeq 30x'].sort_values('sample_id')
        picks = (list(nov[nov.orig_gvcf_dir != ''].sample_id[:2]) + list(nov[nov.orig_gvcf_dir == ''].sample_id[:1]) +
                 list(deep[deep.platform == 'DNBSeq-G400RS 30x'].sort_values('sample_id').sample_id[:1]))
        deep = deep[deep.sample_id.isin(picks)]
        ctrl = ctrl.sort_values('sample_id').head(2)
    else:
        rng = np.random.default_rng(a.seed)
        ctrl = ctrl.iloc[np.sort(rng.choice(len(ctrl), size=min(a.fidelity_controls, len(ctrl)), replace=False))]

    rows = []
    for r in deep.itertuples(index=False):
        rows.append({'sample_id': r.sample_id, 'role': 'deep', 'arm': 'baseline', 'fraction': 1.0, 'seed': 0,
                     'cram': r.cram, 'crai': r.crai, 'platform': r.platform})
        for k in range(1, a.replicates + 1):
            rows.append({'sample_id': r.sample_id, 'role': 'deep', 'arm': f'r{k:02d}', 'fraction': float(r.fraction),
                         'seed': k, 'cram': r.cram, 'crai': r.crai, 'platform': r.platform})
    for r in ctrl.itertuples(index=False):
        rows.append({'sample_id': r.sample_id, 'role': 'fidelity', 'arm': 'baseline', 'fraction': 1.0, 'seed': 0,
                     'cram': r.cram, 'crai': r.crai, 'platform': r.platform})
    plan = pd.DataFrame(rows)
    K.write_tsv(plan, a.out_plan)
    g = pd.concat([deep[deep.orig_gvcf_dir != ''][['sample_id', 'orig_gvcf_dir']],
                   ctrl[['sample_id', 'orig_gvcf_dir']]])
    K.write_tsv(g, a.out_gvcfs)
    K.write_tsv(deep[['sample_id', 'platform', 'batch', 'fraction']], a.out_deep)
    print(f'[arm_plan] {a.mode}: {deep.sample_id.nunique()} deep cases x {a.replicates + 1} arms, '
          f'{len(ctrl)} fidelity controls, {len(g)} original GVCF sets')


if __name__ == '__main__':
    main()
