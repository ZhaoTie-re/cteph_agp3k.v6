#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Omnibus tests of every residue position conditioned on ONE marker,
#           for one case set over the cohorts: each GWAS MHC lead, and each
#           residue signal lead of the case set (stepwise.py, either model).
#           It is the position side of HLA | lead and HLA | signal
#           (conditional.py and hla_given_signal.py cover single markers; the
#           omnibus stepwise conditions positions on whole positions only).
#           The test is omnibus.py's own LRT on the same fixed-model samples,
#           covariates and residue groups, with the marker's dosage (missing
#           calls mean-imputed) as one more covariate: the lead's from the
#           pooled PLINK set, a residue's from the class's PLINK set.
#
#           Check: the scan without the marker is recomputed with the same code
#           and must reproduce OMNIBUS's marginal P of every unit (the script
#           stops otherwise), so the only difference is the marker.
# Inputs  : per cohort (lists in --cohorts order): the class's PLINK set, the
#           fixed-model phenotype/covariates and keep file, the pooled PLINK
#           set; gwas_leads.ids; leads.residue.txt of each model;
#           OMNIBUS's omnibus.<class>.tsv
# Outputs : omnibus_given.<class>.tsv   one row per (cohort, marker, unit)
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import omnibus as O                     # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case-set', required=True)
    p.add_argument('--class', dest='cls', required=True, choices=['residue'])
    p.add_argument('--cohorts', required=True)
    p.add_argument('--bed', nargs='+', required=True)
    p.add_argument('--pheno-fixed', nargs='+', required=True)
    p.add_argument('--keep', nargs='+', required=True)
    p.add_argument('--pooled', nargs='+', required=True, help='pooled .bed per cohort (the GWAS lead)')
    p.add_argument('--leads', required=True, help='gwas_leads.ids')
    p.add_argument('--signal-leads', nargs='*', default=[], help='leads.residue.txt of STEPWISE, either model')
    p.add_argument('--omnibus', required=True, help="OMNIBUS's omnibus.<class>.tsv (the marginal check)")
    p.add_argument('--pheno-col', default='PHENO1')
    p.add_argument('--covar-cols', required=True)
    p.add_argument('--plink2', required=True)
    p.add_argument('--alpha', type=float, default=0.05)
    p.add_argument('--threads', type=int, default=1)
    return p.parse_args()


def lead_dosage(a, pooled, keep, lead, iids, tag):
    """The lead's counted-allele dosage on the fixed-model samples, in iids order; missing calls mean-imputed."""
    prefix = str(Path(pooled).with_suffix(''))
    r = subprocess.run([a.plink2, '--bfile', prefix, '--keep', keep, '--snp', lead, '--export', 'A',
                        '--threads', str(a.threads), '--out', f'lead_{tag}'], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f'ERROR: plink2 --export of {lead} failed for {tag}:\n{r.stdout[-1500:]}')
    raw = pd.read_csv(f'lead_{tag}.raw', sep='\t', dtype={'FID': str, 'IID': str}, na_values=['NA']).set_index('IID')
    if raw.shape[1] != 6 or set(raw.index) != set(iids):
        sys.exit(f'ERROR: {tag}: {lead} missing from the pooled set, or its samples differ from the fixed-model samples')
    g = raw.iloc[:, -1].reindex(iids).to_numpy(dtype=float)
    return np.where(np.isnan(g), np.nanmean(g), g), int(np.isnan(g).sum())


def main():
    a = parse_args()
    cohorts = a.cohorts.split(',')
    if not (len(a.bed) == len(a.pheno_fixed) == len(a.keep) == len(a.pooled) == len(cohorts)):
        sys.exit('ERROR: one PLINK set, phenotype file, keep file and pooled set per cohort')
    leads = [ln.strip() for ln in open(a.leads) if ln.strip()]
    signals = list(dict.fromkeys(ln.strip() for f in a.signal_leads for ln in open(f) if ln.strip()))
    marg = pd.read_csv(a.omnibus, sep='\t', dtype={'unit': str})
    cov = a.covar_cols.split(',')
    rows = []
    for c, bed, ph, kp, pooled in zip(cohorts, a.bed, a.pheno_fixed, a.keep, a.pooled):
        fx, x = O.read_cohort(a, bed, ph, kp, c)
        y = fx[a.pheno_col].astype(float).to_numpy()
        C = np.column_stack([np.ones(len(fx)), fx[cov].astype(float).to_numpy()])
        byu, units = {}, []
        for m in x.columns:
            u = O.unit_of(m, a.cls)
            byu.setdefault(u, []).append(m)
            if u not in units:
                units.append(u)
        cats = {}
        for u, ms in byu.items():
            xu = x[ms].rename(columns=lambda m: O.label_of(m, a.cls))
            cats[u] = O.categories(xu, ~xu.isna().any(axis=1).to_numpy())
        thr = a.alpha / sum(v is not None for v in cats.values())
        # the marginal scan, recomputed, must equal OMNIBUS's
        m0 = O.scan(y, C, cats, [], units)
        mo = marg[marg.cohort == c].set_index('unit')
        bad = [u for u in units if np.isfinite(m0[u]['p']) != np.isfinite(mo.at[u, 'p']) or
               (np.isfinite(m0[u]['p']) and not np.isclose(m0[u]['p'], mo.at[u, 'p'], rtol=1e-4, atol=0))]
        if bad:
            sys.exit(f'ERROR: {c}: the recomputed marginal omnibus P differs from OMNIBUS for {bad[:5]}')
        given = [(ld, 'gwas_lead') for ld in leads] + [(m, 'signal_lead') for m in signals if m in x.columns]
        for mk, why in given:
            if why == 'gwas_lead':
                g, n_miss = lead_dosage(a, pooled, kp, mk, fx.IID.tolist(), c)
            else:
                g = x[mk].to_numpy(dtype=float)
                n_miss = int(np.isnan(g).sum())
                g = np.where(np.isnan(g), np.nanmean(g), g)
            m1 = O.scan(y, np.column_stack([C, g]), cats, [], units, thr)
            for u in units:
                r1, r0 = m1[u], m0[u]
                rows.append({'case_set': a.case_set, 'class': a.cls, 'cohort': c, 'conditioned_on': mk, 'why': why,
                             'unit': u, 'df': r1.get('df', np.nan), 'chisq': r1.get('chisq', np.nan), 'p': r1['p'],
                             'p_marginal': r0['p'], 'threshold': thr, 'sig_bonf': bool(r1.get('sig', False)),
                             'missing_imputed': n_miss, 'errcode': r1['errcode']})
        print(f'[omnibus_given] {a.case_set} {c}: {len(units)} units, marginal reproduced; given '
              f"{', '.join(m for m, _ in given)}")
    pd.DataFrame(rows).to_csv(f'omnibus_given.{a.cls}.tsv', sep='\t', index=False, float_format='%.6g')


if __name__ == '__main__':
    main()
