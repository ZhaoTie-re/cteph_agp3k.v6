#!/usr/bin/env python3
"""Platform subsets of one cohort, with ALL the cohort's controls in every subset (user decision
2026-10-30: shared controls).

  pooled          every case (the published analysis; its files must equal assoc_hla.v2's 00.prep)
  only_<slug>     the cases of one platform (a stratum of the by-platform analysis)
  without_<slug>  every case except one platform's (leave-one-platform-out)

The pooled files come from assoc_hla.v2's own prep_pheno_cov.py (--case-set all); a subset keeps
their rows of the subset's samples, unchanged, so formats and covariates are identical. Written per
subset: the random-model table (SAIGE step 1), the fixed-model table and keep list (unrelated
samples), the rare-variant include list (tuning.rv call-set samples, for rvtest), and the counts.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohort', required=True)
    ap.add_argument('--samples', required=True, help='sample_table.tsv')
    ap.add_argument('--random', required=True, help='pooled pheno_covar.random.tsv (prep_pheno_cov.py)')
    ap.add_argument('--fixed', required=True, help='pooled pheno_covar.fixed.tsv')
    ap.add_argument('--keep', required=True, help='pooled fixed.keep')
    ap.add_argument('--published-prep', required=True, help="assoc_hla.v2 results/all/<cohort>/00.prep")
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # the pooled files must be the published ones (gate G5, first half)
    same = {n: md5(p) == md5(Path(a.published_prep) / n) for n, p in
            (('pheno_covar.random.tsv', a.random), ('pheno_covar.fixed.tsv', a.fixed), ('fixed.keep', a.keep))}
    if not all(same.values()):
        sys.exit(f'ABORT: the pooled prep files differ from the published ones: {same}')

    st = K.read_tsv(a.samples).set_index('sample_id')
    rnd = pd.read_csv(a.random, sep='\t', dtype=str)
    fix = pd.read_csv(a.fixed, sep='\t', dtype=str)
    keep = pd.read_csv(a.keep, sep='\t', dtype=str, header=None)
    rv = set(pd.read_csv(K.TUNING_RV / a.cohort / '02.callset_filter' / 'filtered.fam', sep=r'\s+',
                         header=None, dtype=str)[1])
    plat = st.platform.to_dict()
    is_case = rnd.set_index('IID').PHENO1.eq('1').to_dict()
    if set(is_case) - set(plat):
        sys.exit('ABORT: phenotype samples missing from the sample table')
    if any(plat[s] != K.CONTROL_PLATFORM for s, c in is_case.items() if not c):
        sys.exit('ABORT: a control is not on the control platform')

    subsets = {'pooled': lambda s: True}
    for pf in K.CASE_PLATFORMS:
        slug = K.PLATFORM_SLUG[pf]
        subsets[f'only_{slug}'] = (lambda pf: lambda s: (not is_case[s]) or plat[s] == pf)(pf)
        subsets[f'without_{slug}'] = (lambda pf: lambda s: (not is_case[s]) or plat[s] != pf)(pf)

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
        counts.append({'cohort': a.cohort, 'subset': name,
                       'cases_random': int(r.PHENO1.eq('1').sum()), 'controls_random': int(r.PHENO1.eq('0').sum()),
                       'cases_fixed': int(f.PHENO1.eq('1').sum()), 'controls_fixed': int(f.PHENO1.eq('0').sum()),
                       'cases_rv': int(sum(is_case[s] for s in inc)), 'controls_rv': int(sum(not is_case[s] for s in inc))})
    c = pd.DataFrame(counts)
    # every stratum keeps the same controls; strata cases add up to the pooled cases
    strata = c[c.subset.str.startswith('only_')]
    pooled = c[c.subset == 'pooled'].iloc[0]
    for col in ('cases_random', 'cases_fixed', 'cases_rv'):
        if strata[col].sum() != pooled[col]:
            sys.exit(f'ABORT: strata {col} do not add up to the pooled count')
    for col in ('controls_random', 'controls_fixed', 'controls_rv'):
        if c[col].nunique() != 1:
            sys.exit(f'ABORT: subsets differ in {col}')
    c.to_csv(out / 'subset_counts.tsv', sep='\t', index=False)
    (out / 'subsets.json').write_text(json.dumps({'pooled_equals_published': same, 'subsets': list(subsets)},
                                                 indent=2) + '\n')
    print(c.to_string(index=False))


if __name__ == '__main__':
    main()
