#!/usr/bin/env python3
"""The number of replicates R, from the pilot's measured cost (rule fixed in the plan, 2026-10-30).

Cost is billed core-hours: realtime x the cores the process is allocated (nextflow.config). From the
pilot trace (COMPLETED tasks, or CACHED ones that carry their original realtime):
  per genome     EXTRACT_READS
  per genome-arm SUBSAMPLE, CALL_GVCF, HLA_SUBSAMPLE, HLAHD
  per arm        SWAP, RETEST_GWAS, RETEST_GENE, HLA_INPUTS, HLA_CHAIN, HLA_MATRICES, RETEST_HLA,
                 RETEST_OMNIBUS
  per joint call JOINT, GENOTYPE_QC (scaled by the number of GVCFs they read)
projected for the full design (89 deep cases x (R + 1) arms, fidelity controls at f = 1, the surviving
original GVCFs). R = the largest value in [3, 10] with cost <= the budget (1,000 core-h) and an
estimated wall time <= 48 h at the concurrency given (default 200 cores held at once).
With --retest, the pilot's replicate spread is added: z of the GWAS leads (pooled, best_guess) and of
DQB1:57:V per replicate minus the baseline (only the pilot genomes are down-sampled, so it is small).
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

CORES = {'EXTRACT_READS': 2, 'SUBSAMPLE': 2, 'CALL_GVCF': 2, 'HLA_SUBSAMPLE': 2, 'HLAHD': 4, 'JOINT': 4,
         'GENOTYPE_QC': 4, 'SWAP': 2, 'RETEST_GWAS': 2, 'RETEST_GENE': 2, 'HLA_INPUTS': 2, 'HLA_CHAIN': 8,
         'HLA_MATRICES': 2, 'RETEST_HLA': 2, 'RETEST_OMNIBUS': 2, 'CONCAT_QC': 1, 'CONCORDANCE': 2}
PER_GENOME_ARM = ['SUBSAMPLE', 'CALL_GVCF', 'HLA_SUBSAMPLE', 'HLAHD']
PER_ARM = ['SWAP', 'RETEST_GWAS', 'RETEST_GENE', 'HLA_INPUTS', 'HLA_CHAIN', 'HLA_MATRICES', 'RETEST_HLA',
           'RETEST_OMNIBUS']


def seconds(x):
    if not isinstance(x, str) or x in ('-', ''):
        return np.nan
    t = 0.0
    for v, u in re.findall(r'([\d.]+)\s*(ms|h|m|s|d)', x):
        t += float(v) * {'ms': 1e-3, 's': 1, 'm': 60, 'h': 3600, 'd': 86400}[u]
    return t


def replicate_dz(retest, label, signals):
    """Per cohort x marker (every re-tested GWAS variant, and DQB1:57:V): z at the baseline, and the mean
    and SD over replicates of z - z(baseline); lead = a GWAS lead or DQB1:57:V."""
    s = pd.read_csv(signals, sep='\t', dtype=str)
    leads = set(s[(s.kind == 'gwas') & (s.role == 'Lead')].variant_id) | {'DQB1:57:V'}
    rows = []
    for d in sorted(Path(retest, label).glob('*/*_mainland')):
        arm, cohort = d.parent.name, d.name
        g = pd.read_csv(d / 'gwas.best_guess.tsv', sep='\t', dtype={'variant_id': str})
        g = g[g.subset == 'pooled']
        rows += [{'arm': arm, 'cohort': cohort, 'marker': r.variant_id, 'z': r.z} for r in g.itertuples()]
        h = pd.read_csv(d / 'assoc.random.residue.tsv', sep='\t', dtype={'marker': str})
        h = h[h.marker == 'DQB1:57:V']
        rows += [{'arm': arm, 'cohort': cohort, 'marker': 'DQB1:57:V', 'z': r.beta / r.se} for r in h.itertuples()]
    z = pd.DataFrame(rows).pivot_table(index=['cohort', 'marker'], columns='arm', values='z')
    dz = z[[c for c in z.columns if c != 'baseline']].sub(z['baseline'], axis=0)
    out = pd.DataFrame({'z_baseline': z['baseline'], 'dz_mean': dz.mean(axis=1),
                        'dz_sd': dz.std(axis=1, ddof=1)}).reset_index()
    out['lead'] = out.marker.isin(leads)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--traces', nargs='+', required=True)
    ap.add_argument('--pilot-plan', required=True)
    ap.add_argument('--pilot-gvcfs', required=True, help="the pilot's orig_gvcfs.tsv")
    ap.add_argument('--deep', type=int, default=89)
    ap.add_argument('--fidelity-controls', type=int, default=20)
    ap.add_argument('--original-gvcfs', type=int, default=32)
    ap.add_argument('--chromosomes', type=int, default=14)
    ap.add_argument('--budget', type=float, default=1000.0)
    ap.add_argument('--max-wall-h', type=float, default=48.0)
    ap.add_argument('--cores-held', type=int, default=200)
    ap.add_argument('--retest', help='results/06.retest (optional: the replicate spread of z)')
    ap.add_argument('--label', default='pilot')
    ap.add_argument('--signals', help='results/02.targets/signal_variants.tsv (with --retest)')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    t = pd.concat([pd.read_csv(f, sep='\t', dtype=str) for f in a.traces], ignore_index=True)
    t = t[t.status.isin(['COMPLETED', 'CACHED'])].copy()
    t['process'] = t.name.str.split(' ').str[0]
    t['sec'] = t.realtime.map(seconds)
    t = t.drop_duplicates('hash', keep='last')
    t = t[t.process.isin(CORES)]
    med = t.groupby('process').sec.median()
    n = t.groupby('process').size()
    core_h = {p: med[p] / 3600 * CORES[p] for p in med.index}

    plan = pd.read_csv(a.pilot_plan, sep='\t', dtype=str)
    pilot_gvcfs = len(plan)                                    # one GVCF per genome-arm of the pilot
    orig_pilot = len(pd.read_csv(a.pilot_gvcfs, sep='\t', dtype=str))

    rows = []
    for R in range(3, 11):
        arms = R + 1
        genome_arms = a.deep * arms
        gvcfs = genome_arms + a.fidelity_controls + a.original_gvcfs
        scale = gvcfs / (pilot_gvcfs + orig_pilot)
        c = 0.0
        c += core_h.get('EXTRACT_READS', 0) * (a.deep + a.fidelity_controls)
        c += sum(core_h.get(p, 0) for p in ('SUBSAMPLE', 'CALL_GVCF')) * (genome_arms + a.fidelity_controls)
        c += sum(core_h.get(p, 0) for p in ('HLA_SUBSAMPLE', 'HLAHD')) * genome_arms
        per_arm = (core_h.get('SWAP', 0) * 3 + core_h.get('RETEST_GWAS', 0) * 6 + core_h.get('RETEST_GENE', 0) * 6 +
                   core_h.get('HLA_INPUTS', 0) + core_h.get('HLA_CHAIN', 0) + core_h.get('HLA_MATRICES', 0) * 6 +
                   core_h.get('RETEST_HLA', 0) * 6 + core_h.get('RETEST_OMNIBUS', 0))
        c += per_arm * arms
        joint = (core_h.get('JOINT', 0) + core_h.get('GENOTYPE_QC', 0)) * a.chromosomes * scale
        c += joint
        # wall: the parallel part at the cores held, plus the serial tail (the slowest joint call scaled,
        # then one HLA chain per arm run side by side)
        crit = (med.get('EXTRACT_READS', 0) + med.get('SUBSAMPLE', 0) + med.get('CALL_GVCF', 0) +
                med.get('JOINT', 0) * scale + med.get('GENOTYPE_QC', 0) + med.get('HLA_CHAIN', 0)) / 3600
        wall = c / a.cores_held + crit
        rows.append({'R': R, 'arms': arms, 'gvcfs_in_joint_call': gvcfs, 'core_h': round(c, 1), 'wall_h_est': round(wall, 2)})
    proj = pd.DataFrame(rows)
    ok = proj[(proj.core_h <= a.budget) & (proj.wall_h_est <= a.max_wall_h)]
    R = int(ok.R.max()) if len(ok) else 3
    rep = {'R': R, 'rule': f'largest R in [3, 10] with cost <= {a.budget} core-h and wall <= {a.max_wall_h} h '
                            f'at {a.cores_held} cores held',
           'median_realtime_s': {p: round(float(v), 1) for p, v in med.items()},
           'tasks_measured': {p: int(v) for p, v in n.items()},
           'core_h_per_task': {p: round(float(v), 4) for p, v in core_h.items()},
           'projection': proj.to_dict('records')}
    if a.retest:
        dz = replicate_dz(a.retest, a.label, a.signals)
        moved = dz.dz_sd > 0
        rep['pilot_delta_z'] = {'markers_by_cohort': dz.groupby('cohort').size().to_dict(),
                                'share_moved_by_cohort': moved.groupby(dz.cohort).mean().round(4).to_dict(),
                                'median_sd_moved_by_cohort': dz[moved].groupby('cohort').dz_sd.median().round(4).to_dict(),
                                'max_abs_mean_by_cohort': dz.dz_mean.abs().groupby(dz.cohort).max().round(4).to_dict(),
                                'leads': dz[dz.lead].round(4).to_dict('records')}
        print(json.dumps({k: v for k, v in rep['pilot_delta_z'].items() if k != 'leads'}, indent=1))
    Path(a.out).write_text(json.dumps(rep, indent=2) + '\n')
    print(proj.to_string(index=False))
    print(f'R = {R}')


if __name__ == '__main__':
    main()
