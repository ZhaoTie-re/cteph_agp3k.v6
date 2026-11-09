#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : MAIN ANALYSIS. Omnibus tests of the residue positions of one case
#           set in the FIXED model, the three cohorts jointly: each amino-acid
#           position is tested as a whole, and independent positions are
#           selected stepwise, conditioning on all categories of every entered
#           position. Residues only: gene-level tests over the alleles of a gene
#           are not run (not standard in HLA fine-mapping, and their 'unlisted'
#           category differs between cases and controls in every gene, as a
#           platform difference in typing would; METHODS §8). The allele
#           branches below are kept only because units are named by class.
#
#           UNIT AND CATEGORIES. A residue position's categories are its
#           residues; a gene's, its 2-field alleles. Copies carrying a residue or
#           allele the matrix does not list (2 minus the listed copies) form one
#           more category, 'unlisted'. Within the cohort's fixed-model samples
#           called at the unit: a category with copy frequency < 1 % joins the
#           pool 'rare'; a pool still < 0.5 % (the single-marker floor) joins the
#           reference; the reference is the most frequent category. df =
#           categories - 1. A unit with no category left besides the reference
#           is not tested.
#
#           TEST. Logistic regression on the case set's fixed-model (unrelated)
#           samples called at the unit and at every conditioning unit:
#           PHENO1 ~ 1 + SEX + PC1-10 [+ the entered units' categories] + the
#           unit's categories, against the same model without the unit's:
#           likelihood-ratio test, chi-square with df = the rank the unit's
#           columns add (a column collinear with the model adds none). Copies
#           are additive (0/1/2), as in the single-marker scans.
#
#           THRESHOLD. Bonferroni: 0.05 / the units tested in the cohort's
#           marginal scan. A fit that separated (|log OR| > 15) or did not
#           converge keeps its P but is never significant.
#
#           STEPWISE, the rule of STEPWISE for single markers: candidates = the
#           units significant in at least one cohort; the lead is the candidate
#           whose largest P over the cohorts is smallest (a unit a cohort did not
#           test counts P = 1 there); the walk stops when no unit is a candidate,
#           after at most --max-signals signals. A signal's MEMBERS are the units
#           significant before its lead entered and in no cohort after, the lead
#           included. CONSISTENCY is the number of cohorts in which a unit is
#           significant (an omnibus test has no single direction).
# Inputs  : per cohort (lists in --cohorts order): the class's PLINK set
#           (BUILD_MARKERS), pheno_covar.fixed.tsv, fixed.keep (PREP)
# Outputs : omnibus.<class>.tsv          marginal: every unit x cohort, with the
#                                        categories and their ORs
#           omnibus_stepwise.<class>.tsv every unit x cohort x step
#           omnibus_steps.<class>.tsv    one row per step
#           omnibus_signals.<class>.tsv  one row per entered unit
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import vocab as V                      # noqa: E402

FREQ_OWN = 0.01                        # a category below this joins 'rare'
FREQ_POOL = 0.005                      # a 'rare' pool below this joins the reference


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case-set', required=True)
    p.add_argument('--class', dest='cls', required=True, choices=['residue'])
    p.add_argument('--cohorts', required=True)
    p.add_argument('--bed', nargs='+', required=True)
    p.add_argument('--pheno-fixed', nargs='+', required=True)
    p.add_argument('--keep', nargs='+', required=True)
    p.add_argument('--pheno-col', default='PHENO1')
    p.add_argument('--covar-cols', required=True)
    p.add_argument('--plink2', required=True)
    p.add_argument('--alpha', type=float, default=0.05)
    p.add_argument('--max-signals', type=int, default=5)
    p.add_argument('--threads', type=int, default=1)
    return p.parse_args()


def unit_of(marker, cls):
    gene, pos, _ = V.parse_marker(marker, cls)
    return f'{gene}:{pos}' if cls == 'residue' else gene


def label_of(marker, cls):
    return V.parse_marker(marker, cls)[2] if cls == 'residue' else marker


def read_cohort(a, bed, pheno, keep, tag):
    """(phenotype/covariate frame, dosage frame of copies present) on the fixed-model samples."""
    fx = pd.read_csv(pheno, sep='\t', dtype={'#FID': str, 'IID': str})
    kp = pd.read_csv(keep, sep=r'\s+', header=None, dtype=str)
    if set(kp[1]) != set(fx.IID) or len(kp) != len(fx):
        sys.exit(f'ERROR: {tag}: fixed.keep and pheno_covar.fixed.tsv differ')
    cols = a.covar_cols.split(',')
    if fx[[a.pheno_col] + cols].isna().any().any() or not set(fx[a.pheno_col].unique()) <= {0, 1}:
        sys.exit(f'ERROR: {tag}: missing or non-0/1 phenotype, or missing covariate')
    prefix = str(Path(bed).with_suffix(''))
    r = subprocess.run([a.plink2, '--bfile', prefix, '--keep', keep, '--export', 'A', '--threads', str(a.threads),
                        '--out', f'dosage_{tag}'], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f'ERROR: plink2 --export failed for {tag}:\n{r.stdout[-1500:]}')
    raw = pd.read_csv(f'dosage_{tag}.raw', sep='\t', dtype={'FID': str, 'IID': str}, na_values=['NA']).set_index('IID')
    bim = pd.read_csv(f'{prefix}.bim', sep='\t', header=None, dtype={1: str})
    x = {}
    for c in raw.columns[5:]:
        name, allele = c.rsplit('_', 1)
        if allele not in ('A', 'P'):
            sys.exit(f'ERROR: plink2 counted allele {allele!r} for {name}')
        v = raw[c].to_numpy(dtype=float)
        x[name] = 2.0 - v if allele == 'A' else v
    if list(x) != bim[1].tolist() or set(raw.index) != set(fx.IID):
        sys.exit(f'ERROR: {tag}: exported markers or samples differ from the .bim or the fixed-model samples')
    return fx, pd.DataFrame(x, index=raw.index).loc[fx.IID]


def categories(xu, called):
    """The unit's categories as copy columns over all rows (0 where the unit is not called), the reference and
    the pooled categories; frequencies over the called rows. None if nothing is left to test."""
    src = xu.fillna(0.0)
    cat = {c: np.where(called, src[c].to_numpy(), 0.0) for c in src.columns}
    unl = np.where(called, 2.0 - src.sum(axis=1).to_numpy(), 0.0)
    if (unl < -1e-9).any():
        sys.exit(f'ERROR: more than 2 copies at the unit of {list(xu.columns)[:3]}')
    if unl.max() > 1e-9:
        cat['unlisted'] = np.clip(unl, 0.0, 2.0)
    n2 = 2.0 * called.sum()
    freq = {k: v.sum() / n2 for k, v in cat.items()}
    own = {k: v for k, v in cat.items() if freq[k] >= FREQ_OWN}
    pool = [k for k in cat if 0 < freq[k] < FREQ_OWN]
    if pool and sum(freq[k] for k in pool) >= FREQ_POOL:
        own['rare'] = sum(cat[k] for k in pool)
    else:
        pool = []
    if len(own) < 2:
        return None
    fq = {k: v.sum() / n2 for k, v in own.items()}
    ref = max(fq, key=fq.get)
    return {'test': {k: v for k, v in own.items() if k != ref}, 'ref': ref, 'pool': pool, 'freq': fq,
            'called': called}


def logistic(y, X, maxiter=100, tol=1e-9):
    """Logistic MLE by Newton-Raphson with step halving: (coefficients, deviance, converged)."""
    b = np.zeros(X.shape[1])
    dev = np.inf
    for _ in range(maxiter):
        eta = X @ b
        mu = 1 / (1 + np.exp(-eta))
        w = np.clip(mu * (1 - mu), 1e-12, None)
        step = np.linalg.lstsq(X.T @ (X * w[:, None]), X.T @ (y - mu), rcond=None)[0]
        t = 1.0
        while True:
            nb = b + t * step
            ne = np.clip(X @ nb, -50, 50)
            nd = -2 * np.sum(y * ne - np.logaddexp(0, ne))
            if nd <= dev + 1e-12 or t < 1e-6:
                break
            t /= 2
        b, conv = nb, abs(dev - nd) < tol * (1 + abs(nd))
        dev = nd
        if conv:
            return b, dev, True
    return b, dev, False


def independent(X0, Xa):
    """Columns of Xa that add rank to X0 (None: nothing), in order."""
    keep, cur = [], X0
    r = 0 if cur is None else np.linalg.matrix_rank(cur)
    for j in range(Xa.shape[1]):
        nxt = Xa[:, [j]] if cur is None else np.column_stack([cur, Xa[:, j]])
        rn = np.linalg.matrix_rank(nxt)
        if rn > r:
            keep.append(j)
            cur, r = nxt, rn
    return keep


def test_unit(y, C, unit_cats, cond_cats, ok):
    """LRT of one unit given covariates C and the conditioning units, on the rows ok."""
    names = list(unit_cats['test'])
    Xa = np.column_stack([unit_cats['test'][k] for k in names])
    X0 = C
    for cc in cond_cats:
        X0 = np.column_stack([X0] + [v for v in cc['test'].values()])
    X0, Xa, yy = X0[ok], Xa[ok], y[ok]
    cols = independent(X0, Xa)
    out = {'n': int(ok.sum()), 'cases': int(yy.sum()), 'controls': int(len(yy) - yy.sum()), 'df': len(cols)}
    if not cols:
        return {**out, 'chisq': np.nan, 'p': np.nan, 'errcode': 'COLLINEAR', 'ors': ''}
    keep0 = independent(None, X0)
    X0 = X0[:, keep0]
    _, d0, c0 = logistic(yy, X0)
    b1, d1, c1 = logistic(yy, np.column_stack([X0, Xa[:, cols]]))
    chisq = max(d0 - d1, 0.0)
    beta = b1[X0.shape[1]:]
    err = '.' if (c0 and c1) else 'UNFINISHED'
    if np.abs(beta).max() > 15:
        err = 'SEPARATION' if err == '.' else err
    ors = ','.join(f'{names[j]}:{np.exp(bb):.4g}' for j, bb in zip(cols, beta))
    return {**out, 'chisq': chisq, 'p': float(stats.chi2.sf(chisq, len(cols))), 'errcode': err, 'ors': ors}


def scan(y, C, cats, entered, units, sig_thr=None):
    """Every unit given the entered ones, for one cohort; rows keyed by unit."""
    rows = {}
    cond = [cats[e] for e in entered if cats.get(e) is not None]
    for u in units:
        if u in entered:
            rows[u] = {'errcode': 'ENTERED', 'p': np.nan}
            continue
        uc = cats.get(u)
        if uc is None:
            rows[u] = {'errcode': 'NOT_TESTED', 'p': np.nan}
            continue
        ok = uc['called'].copy()
        for cc in cond:
            ok &= cc['called']
        r = test_unit(y, C, uc, cond, ok)
        if sig_thr is not None:                               # a fit that separated or did not finish never counts
            r['sig'] = bool(np.isfinite(r['p']) and r['p'] < sig_thr and r['errcode'] == '.')
        rows[u] = r
    return rows


def main():
    a = parse_args()
    cohorts = a.cohorts.split(',')
    if not (len(a.bed) == len(a.pheno_fixed) == len(a.keep) == len(cohorts)):
        sys.exit('ERROR: one PLINK set, phenotype file and keep file per cohort')
    cov = a.covar_cols.split(',')
    data, units = {}, []
    for c, bed, ph, kp in zip(cohorts, a.bed, a.pheno_fixed, a.keep):
        fx, x = read_cohort(a, bed, ph, kp, c)
        y = fx[a.pheno_col].astype(float).to_numpy()
        C = np.column_stack([np.ones(len(fx)), fx[cov].astype(float).to_numpy()])
        byu = {}
        for m in x.columns:
            byu.setdefault(unit_of(m, a.cls), []).append(m)
        cats = {}
        for u, ms in byu.items():
            xu = x[ms].rename(columns=lambda m: label_of(m, a.cls))
            blank = xu.isna()
            if (blank.any(axis=1) & ~blank.all(axis=1)).any():
                sys.exit(f'ERROR: {c} {u}: some samples are blank in only part of the unit')
            cats[u] = categories(xu, ~blank.any(axis=1).to_numpy())
            if u not in units:
                units.append(u)
        data[c] = {'y': y, 'C': C, 'cats': cats, 'tested_units': sum(v is not None for v in cats.values())}

    # ── step 0: the marginal scans ──────────────────────────────────────────
    thr = {c: a.alpha / d['tested_units'] for c, d in data.items()}
    marg = {c: scan(d['y'], d['C'], d['cats'], [], units, thr[c]) for c, d in data.items()}
    rows = []
    for c in cohorts:
        for u in units:
            r, k = marg[c][u], data[c]['cats'].get(u)
            gene = u.split(':')[0]
            rows.append({'case_set': a.case_set, 'class': a.cls, 'cohort': c, 'unit': u, 'gene': gene,
                         'position': u.split(':')[1] if a.cls == 'residue' else '',
                         'reference': k['ref'] if k else '', 'categories': ','.join(k['test']) if k else '',
                         'pooled': ','.join(k['pool']) if k else '', 'n': r.get('n', np.nan),
                         'cases': r.get('cases', np.nan), 'controls': r.get('controls', np.nan),
                         'df': r.get('df', np.nan), 'chisq': r.get('chisq', np.nan), 'p': r['p'],
                         'threshold': thr[c], 'sig_bonf': bool(r.get('sig', False)), 'ors': r.get('ors', ''),
                         'errcode': r['errcode'] if k else 'NOT_TESTED'})
    pd.DataFrame(rows).to_csv(f'omnibus.{a.cls}.tsv', sep='\t', index=False, float_format='%.6g')

    # ── stepwise ────────────────────────────────────────────────────────────
    long, steps, signals = [], [], []
    entered, cur, step = [], marg, 0

    def record(step, res):
        for c in cohorts:
            for u in units:
                r = res[c][u]
                long.append({'case_set': a.case_set, 'class': a.cls, 'step': step,
                             'conditioned_on': ';'.join(entered) or '-', 'cohort': c, 'unit': u,
                             'df': r.get('df', np.nan), 'chisq': r.get('chisq', np.nan), 'p': r['p'],
                             'sig_bonf': bool(r.get('sig', False)), 'errcode': r['errcode']})

    def sig_any(res, u):
        return [c for c in cohorts if res[c][u].get('sig', False)]

    record(0, cur)
    while True:
        cand = [u for u in units if u not in entered and sig_any(cur, u)]
        row = {'case_set': a.case_set, 'class': a.cls, 'step': step, 'conditioned_on': ';'.join(entered) or '-'}
        if not cand:
            steps.append({**row, 'lead': '', 'entered': False, 'reason': 'no unit significant in any cohort'})
            break
        maxp = {u: max((cur[c][u]['p'] if np.isfinite(cur[c][u]['p']) else 1.0) for c in cohorts) for u in cand}
        lead = min(cand, key=lambda u: (maxp[u], min(cur[c][u]['p'] for c in cohorts if np.isfinite(cur[c][u]['p']))))
        steps.append({**row, 'lead': lead, 'entered': True, 'reason': 'entered', 'lead_max_p': maxp[lead],
                      **{f'p_{c}': cur[c][lead]['p'] for c in cohorts},
                      **{f'n_sig_{c}': sum(cur[c][u].get('sig', False) for u in units) for c in cohorts},
                      'consistency': f'{len(sig_any(cur, lead))}/{len(cohorts)}'})
        before = {u for u in units if sig_any(cur, u)}
        entered.append(lead)
        step += 1
        nxt = {c: scan(d['y'], d['C'], d['cats'], entered, units, thr[c]) for c, d in data.items()}
        record(step, nxt)
        after = {u for u in units if u not in entered and sig_any(nxt, u)}
        members = [u for u in units if u in before and u not in after]
        signals.append({'case_set': a.case_set, 'class': a.cls, 'signal': len(entered), 'lead': lead,
                        'step': step - 1, 'conditioned_on': ';'.join(entered[:-1]) or '-',
                        'consistency': f'{len(sig_any(cur, lead))}/{len(cohorts)}',
                        'consistency_marginal': f'{len(sig_any(marg, lead))}/{len(cohorts)}',
                        **{f'p_{c}': cur[c][lead]['p'] for c in cohorts},
                        **{f'p_marginal_{c}': marg[c][lead]['p'] for c in cohorts},
                        **{f'df_{c}': cur[c][lead].get('df', np.nan) for c in cohorts},
                        'n_members': len(members), 'members': ';'.join(members)})
        cur = nxt
        if len(entered) >= a.max_signals:
            steps.append({'case_set': a.case_set, 'class': a.cls, 'step': step, 'conditioned_on': ';'.join(entered),
                          'lead': '', 'entered': False, 'reason': f'{a.max_signals} signals entered'})
            break
    pd.DataFrame(long).to_csv(f'omnibus_stepwise.{a.cls}.tsv', sep='\t', index=False, float_format='%.6g')
    pd.DataFrame(steps).to_csv(f'omnibus_steps.{a.cls}.tsv', sep='\t', index=False, float_format='%.6g')
    sig_cols = ['case_set', 'class', 'signal', 'lead', 'step', 'conditioned_on', 'consistency', 'consistency_marginal'] + \
        [f'{k}_{c}' for k in ('p', 'p_marginal', 'df') for c in cohorts] + ['n_members', 'members']
    pd.DataFrame(signals, columns=sig_cols).to_csv(f'omnibus_signals.{a.cls}.tsv', sep='\t', index=False,
                                                   float_format='%.6g')
    print(f"[omnibus] {a.case_set} {a.cls}: units tested " +
          ', '.join(f"{c} {data[c]['tested_units']} (Bonferroni {thr[c]:.3g})" for c in cohorts) +
          f"; signals: {', '.join(s['lead'] + ' ' + s['consistency'] for s in signals) or 'none'}")


if __name__ == '__main__':
    main()
