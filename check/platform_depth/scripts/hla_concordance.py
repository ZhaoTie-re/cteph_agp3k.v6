#!/usr/bin/env python3
"""HLA concordance of the re-typed deep cases (user request 2026-10-30: a reported result).

From hla.typing.v2's own call table (02.evidence/call_evidence.tsv) of the archive (= the typing the
association used) and of each arm's side run:

  baseline_vs_archive    f = 1 re-typing of the archived reads vs the archive (gate G9: HLA-HD and the
                         snapshot reproduce the typing)
  replicate_vs_baseline  each replicate vs the baseline (the depth effect)
  replicate_vs_archive   each replicate vs the archive (the total change)

Per sample x gene: the call state of each (heterozygous / homozygous / present / absent / ambiguous /
unresolved / failed / missing), the two-field allele copies shared by the two genotypes (0-2, unordered)
when both are typed, and the allele read support. Summaries per gene and arm. Also the dosage of the
ORIGINAL tested allele and residue columns, arm against archive, for the deep cases (DQB1:57:V itemised),
and the read support of the deep cases after down-sampling against the controls' (the HLA depth check).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TYPED = {'heterozygous', 'homozygous', 'present', 'absent'}


def evidence(path, ids):
    e = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    return e[e.sample_id.isin(ids)].set_index(['sample_id', 'gene'])


def shared_copies(a, b):
    x = sorted(v for v in a if v)
    y = sorted(v for v in b if v)
    n = 0
    y = list(y)
    for v in x:
        if v in y:
            y.remove(v)
            n += 1
    return n


def compare(ref, test, comp, arm):
    rows = []
    for key in ref.index.intersection(test.index):
        r, t = ref.loc[key], test.loc[key]
        both = r.state in TYPED and t.state in TYPED
        ra, ta = (r.allele_1, r.allele_2), (t.allele_1, t.allele_2)
        n_cop = max(len([v for v in ra if v]), len([v for v in ta if v]))
        rows.append({'comparison': comp, 'arm': arm, 'sample_id': key[0], 'gene': key[1],
                     'state_ref': r.state, 'state_test': t.state, 'alleles_ref': '/'.join(ra),
                     'alleles_test': '/'.join(ta), 'both_typed': int(both),
                     'copies': n_cop if both else 0, 'copies_shared': shared_copies(ra, ta) if both else 0,
                     'identical': int(r.state == t.state and sorted(ra) == sorted(ta)),
                     'reads_ref': pd.to_numeric(r.reads_1, errors='coerce') + pd.to_numeric(r.reads_2, errors='coerce'),
                     'reads_test': pd.to_numeric(t.reads_1, errors='coerce') + pd.to_numeric(t.reads_2, errors='coerce'),
                     'confidence_ref': r.confidence, 'confidence_test': t.confidence})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--archive', required=True, help='hla.typing.v2 02.evidence/call_evidence.tsv')
    ap.add_argument('--arms', nargs='+', required=True, help='<arm>=<side 02.evidence/call_evidence.tsv>')
    ap.add_argument('--deep', required=True, help='deep_cases.tsv')
    ap.add_argument('--samples', required=True, help='sample_table.tsv')
    ap.add_argument('--matrices', nargs='*', default=[], help='<arm>:<cohort>:<class>=<dosage.tsv>')
    ap.add_argument('--original-matrices', required=True, help='hla.typing.v2 06.assoc_inputs directory')
    ap.add_argument('--out-prefix', required=True)
    a = ap.parse_args()

    deep = pd.read_csv(a.deep, sep='\t', dtype=str)
    ids = set(deep.sample_id)
    arch = evidence(a.archive, ids)
    arms = dict(x.split('=', 1) for x in a.arms)
    ev = {arm: evidence(p, ids) for arm, p in arms.items()}
    rows = []
    if 'baseline' in ev:
        rows += compare(arch, ev['baseline'], 'baseline_vs_archive', 'baseline')
    for arm in sorted(k for k in ev if k != 'baseline'):
        if 'baseline' in ev:
            rows += compare(ev['baseline'], ev[arm], 'replicate_vs_baseline', arm)
        rows += compare(arch, ev[arm], 'replicate_vs_archive', arm)
    long = pd.DataFrame(rows)
    long.to_csv(f'{a.out_prefix}.calls_long.tsv', sep='\t', index=False)

    g = long.groupby(['comparison', 'arm', 'gene'])
    summ = pd.DataFrame({
        'samples': g.sample_id.nunique(), 'identical_share': g.identical.mean(),
        'both_typed': g.both_typed.sum(), 'copies': g.copies.sum(), 'copies_shared': g.copies_shared.sum(),
        'state_changed': g.apply(lambda d: int((d.state_ref != d.state_test).sum()), include_groups=False),
        'reads_ref_median': g.reads_ref.median(), 'reads_test_median': g.reads_test.median()}).reset_index()
    summ['allele_concordance'] = summ.copies_shared / summ.copies.where(summ.copies > 0)
    tot = long.groupby(['comparison', 'arm']).agg(samples=('sample_id', 'nunique'), identical_share=('identical', 'mean'),
                                                  copies=('copies', 'sum'), copies_shared=('copies_shared', 'sum')).reset_index()
    tot['gene'] = 'all'
    tot['allele_concordance'] = tot.copies_shared / tot.copies
    summ = pd.concat([summ, tot], ignore_index=True)
    summ.to_csv(f'{a.out_prefix}.summary.tsv', sep='\t', index=False)

    # read support: deep cases per arm vs the controls of the archive (full mainland)
    st = pd.read_csv(a.samples, sep='\t', dtype=str)
    ctrl = set(st[(st.status == 'control') & (st.in_full_mainland == '1')].sample_id)
    allarch = pd.read_csv(a.archive, sep='\t', dtype=str, keep_default_na=False)
    allarch['reads'] = pd.to_numeric(allarch.reads_1, errors='coerce') + pd.to_numeric(allarch.reads_2, errors='coerce')
    cr = allarch[allarch.sample_id.isin(ctrl)].groupby('gene').reads.median().rename('controls_median')
    sup = [arch.reset_index().assign(arm='archive')] + [e.reset_index().assign(arm=k) for k, e in ev.items()]
    sup = pd.concat(sup)
    sup['reads'] = pd.to_numeric(sup.reads_1, errors='coerce') + pd.to_numeric(sup.reads_2, errors='coerce')
    rs = sup.groupby(['arm', 'gene']).reads.median().rename('deep_median').reset_index().merge(cr.reset_index(), on='gene')
    rs['deep_over_controls'] = rs.deep_median / rs.controls_median
    rs.to_csv(f'{a.out_prefix}.read_support.tsv', sep='\t', index=False)

    # dosage of the ORIGINAL tested columns, arm vs archive, deep cases only
    drows = []
    for spec in a.matrices:
        key, path = spec.split('=', 1)
        arm, cohort, cls = key.split(':')
        m = pd.read_csv(path, sep='\t', dtype={'sample_id': str}).set_index('sample_id')
        o = pd.read_csv(Path(a.original_matrices) / cohort / f'{cls}_dosage.tsv', sep='\t',
                        dtype={'sample_id': str}).set_index('sample_id')
        if list(m.columns) != list(o.columns):
            sys.exit(f'ABORT: {key} columns differ from the original matrix')
        sel = [s for s in m.index if s in ids]
        x, y = o.loc[sel], m.loc[sel]
        same = (x.fillna(-1) == y.fillna(-1))
        drows.append({'arm': arm, 'cohort': cohort, 'class': cls, 'samples': len(sel), 'columns': x.shape[1],
                      'cells': int(same.size), 'cells_changed': int((~same).sum().sum()),
                      'to_blank': int((x.notna() & y.isna()).sum().sum()), 'from_blank': int((x.isna() & y.notna()).sum().sum()),
                      'columns_changed': int((~same).any().sum())})
        if cls == 'residue' and 'DQB1:57:V' in m.columns:
            drows[-1]['dqb1_57_v_copies_original'] = int(x['DQB1:57:V'].sum())
            drows[-1]['dqb1_57_v_copies_arm'] = int(y['DQB1:57:V'].sum())
            drows[-1]['dqb1_57_v_samples_changed'] = int((x['DQB1:57:V'].fillna(-1) != y['DQB1:57:V'].fillna(-1)).sum())
    pd.DataFrame(drows).to_csv(f'{a.out_prefix}.dosage.tsv', sep='\t', index=False)
    print(summ[summ.gene == 'all'].round(4).to_string(index=False))


if __name__ == '__main__':
    main()
