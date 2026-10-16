#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Do the ancestry PCs differ between cases and controls? Per cohort,
#           logistic regression of case status on SEX alone against SEX + PC1-K,
#           on the samples of the mixed-model scan: likelihood-ratio chi2(K) and
#           its P, the Wald P of every PC, and McFadden's pseudo-R2 gain.
#           A PC associated with case status is a potential confounder of every
#           variant whose frequency follows that PC.
# Inputs  : assoc_saige <cohort>/00.prep/pheno_covar.tsv (IID PHENO1 SEX PC1_AVG ...)
# Outputs : pc_case_association.tsv
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import pandas as pd
import statsmodels.api as sm
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import csv_list, log  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohorts', required=True)
    ap.add_argument('--pheno', required=True, help='pheno_covar.tsv files, comma-separated, in cohort order')
    ap.add_argument('--n-pcs', type=int, default=10)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    cohorts, files = csv_list(a.cohorts), csv_list(a.pheno)
    assert len(cohorts) == len(files)
    pcs = [f'PC{i}_AVG' for i in range(1, a.n_pcs + 1)]
    rows = []
    for c, f in zip(cohorts, files):
        d = pd.read_csv(f, sep='\t')
        y = d.PHENO1.astype(int)
        assert set(y.unique()) <= {0, 1}
        m0 = sm.Logit(y, sm.add_constant(d[['SEX']].astype(float))).fit(disp=0)
        m1 = sm.Logit(y, sm.add_constant(d[['SEX'] + pcs].astype(float))).fit(disp=0)
        lr = 2 * (m1.llf - m0.llf)
        r = dict(cohort=c, n=len(d), n_case=int(y.sum()), n_ctrl=int((1 - y).sum()), n_pcs=a.n_pcs,
                 lr_chi2=lr, lr_df=a.n_pcs, lr_p=stats.chi2.sf(lr, a.n_pcs),
                 pseudo_r2_sex=m0.prsquared, pseudo_r2_sex_pcs=m1.prsquared,
                 n_pcs_p_lt_0_05=int((m1.pvalues[pcs] < 0.05).sum()))
        r.update({f'p_{p.replace("_AVG", "")}': m1.pvalues[p] for p in pcs})
        rows.append(r)
        log(f'{c}: LR chi2({a.n_pcs}) = {lr:.1f}, P = {r["lr_p"]:.2e}')
    pd.DataFrame(rows).to_csv(a.out, sep='\t', index=False, float_format='%.6g')


if __name__ == '__main__':
    main()
