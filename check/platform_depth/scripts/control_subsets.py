#!/usr/bin/env python3
"""Control-source subsets of one cohort (user request 2026-10-30): every case, with

  ctrlonly_<source>     the controls of one source only
  ctrlwithout_<source>  the controls of every other source (leave one control source out)

Sources are the info sheet's Project column (CONTROL_SOURCES below: Nagahama, BBJ, HTLV-1, ACC, lung
cancer; kept here, not in common.py, which is a cached input of most processes). The rows are those of the
pooled platform subset (make_subsets.py, itself equal to assoc_hla.v2's
00.prep), kept unchanged, so formats and covariates are identical. Written per subset as in make_subsets.py:
random-model table, fixed-model table and keep list, rare-variant include list, and subset_counts.tsv.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402

# The info sheet's Project column, verbatim -> the slug used in subset names
CONTROL_SOURCES = {
    'Control_Nagahama': 'nagahama',
    'Control_BBJ': 'bbj',
    'HTLV-1': 'htlv1',
    'Control_ACC': 'acc',
    'Lung_cancer': 'lungca',
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohort', required=True)
    ap.add_argument('--subsets', required=True, help="the platform subsets folder (pooled.*)")
    ap.add_argument('--info', required=True, help='the info sheet (xlsx)')
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    src, out = Path(a.subsets), Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    rnd = pd.read_csv(src / 'pooled.random.tsv', sep='\t', dtype=str)
    fix = pd.read_csv(src / 'pooled.fixed.tsv', sep='\t', dtype=str)
    keep = pd.read_csv(src / 'pooled.keep', sep='\t', dtype=str, header=None)
    rv = set((src / 'pooled.rv_include.txt').read_text().split())
    info = pd.read_excel(a.info, sheet_name=0, dtype=str)
    info = info.dropna(subset=['ID_JHRPv6'])[['ID_JHRPv6', 'Project']]

    is_case = rnd.set_index('IID').PHENO1.eq('1').to_dict()
    ctrl = [s for s, c in is_case.items() if not c]
    info = info[info.ID_JHRPv6.isin(ctrl)]
    if info.groupby('ID_JHRPv6').Project.nunique().gt(1).any():
        sys.exit('ABORT: a control has two Project values in the info sheet')
    project = info.drop_duplicates('ID_JHRPv6').set_index('ID_JHRPv6').reindex(ctrl).Project
    bad = project[~project.isin(CONTROL_SOURCES)]
    if len(bad):
        sys.exit(f'ABORT: {len(bad)} controls without a known source: {bad.value_counts(dropna=False).to_dict()}')
    source = {s: CONTROL_SOURCES[p] for s, p in project.items()}

    subsets = {}
    for slug in CONTROL_SOURCES.values():
        subsets[f'ctrlonly_{slug}'] = (lambda g: lambda s: is_case[s] or source[s] == g)(slug)
        subsets[f'ctrlwithout_{slug}'] = (lambda g: lambda s: is_case[s] or source[s] != g)(slug)

    counts = []
    for name, keep_fn in subsets.items():
        r = rnd[rnd.IID.map(keep_fn)]
        f = fix[fix.IID.map(keep_fn)]
        k = keep[keep[1].map(keep_fn)]
        r.to_csv(out / f'{name}.random.tsv', sep='\t', index=False)
        f.to_csv(out / f'{name}.fixed.tsv', sep='\t', index=False)
        k.to_csv(out / f'{name}.keep', sep='\t', index=False, header=False)
        inc = [s for s in r.IID if s in rv]
        Path(out / f'{name}.rv_include.txt').write_text('\n'.join(inc) + '\n')
        counts.append({'cohort': a.cohort, 'subset': name, 'source': name.split('_', 1)[1],
                       'cases_random': int(r.PHENO1.eq('1').sum()), 'controls_random': int(r.PHENO1.eq('0').sum()),
                       'cases_fixed': int(f.PHENO1.eq('1').sum()), 'controls_fixed': int(f.PHENO1.eq('0').sum()),
                       'cases_rv': sum(is_case[s] for s in inc), 'controls_rv': sum(not is_case[s] for s in inc)})
    c = pd.DataFrame(counts)
    n_ctrl = len(ctrl)
    only = c[c.subset.str.startswith('ctrlonly_')].set_index('source').controls_random
    wo = c[c.subset.str.startswith('ctrlwithout_')].set_index('source').controls_random
    if only.sum() != n_ctrl or not (wo == n_ctrl - only).all() or c.cases_random.nunique() != 1:
        sys.exit('ABORT: control-source subsets do not add up to the pooled controls')
    c.to_csv(out / 'subset_counts.tsv', sep='\t', index=False)
    print(c.to_string(index=False))


if __name__ == '__main__':
    main()
