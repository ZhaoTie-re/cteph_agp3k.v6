#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : BENCHMARK for one (case set, cohort, class): the exact family-wise
#           threshold by parametric bootstrap on the FIXED-model samples, and
#           the family-wise error rate the GEC and Bonferroni thresholds
#           actually give. No call depends on it.
#
#           NULL PHENOTYPES. The null logistic model PHENO1 ~ 1 + SEX + PC1-10
#           is fitted on the case set's unrelated (fixed-model) samples; each of
#           B replicates draws y* ~ Bernoulli(mu-hat) per sample, so the
#           covariates keep their effect and the genotypes have none. Case
#           labels are not permuted: the PCs predict case status, and a
#           permutation would remove that dependence.
#
#           TWO STATISTICS, on the same replicates:
#             wald   the fixed model's own test: plink2 reruns every test with
#                    ASSOC_FIXED's settings (--glm omit-ref no-firth, the same
#                    covariates, --1, complete-case per marker), in batches of
#                    phenotype columns. The observed PHENO1 is rerun as one more
#                    column, so verify.sh can check the rerun equals ASSOC_FIXED.
#             score  the random model's statistic without the GRM: the logistic
#                    score test on the null model refitted per replicate, blanks
#                    mean-imputed, P from the normal approximation (no SPA), as
#                    SAIGE's p.value.NA.
#           One marker per mirror pair.
#
#           SUMMARY, per statistic s. For replicate b, minP_b = the smallest P
#           over the tests.
#             p_maxt_<s>          alpha quantile of minP: the exact threshold;
#                                 95% interval from 1,000 resamples of minP
#             fwer_<t>_<s>        share of replicates with minP < p_<t> (as the
#                                 calls compare), Wilson 95% interval;
#                                 t = gec, bonf
#           The multivariate-normal max-T of EFFECTIVE_TESTS is copied beside
#           them (p_maxt_mvn, fwer_<t>_mvn) for comparison.
# Inputs  : the class's PLINK set (BUILD_MARKERS), pheno_covar.fixed.tsv and
#           fixed.keep (PREP), tests.<class>.tsv and effective_tests.<class>.tsv
#           (EFFECTIVE_TESTS)
# Outputs : maxt_boot.<class>.tsv           one row
#           maxt_boot.<class>.minp.tsv.gz   one row per replicate
#           maxt_boot.<class>.observed.tsv  the observed phenotype's P per test
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

STATS = ('wald', 'score')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case-set', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--class', dest='cls', required=True, choices=V.CLASSES)
    p.add_argument('--bfile', required=True, help="the class's PLINK prefix (BUILD_MARKERS)")
    p.add_argument('--pheno-fixed', required=True, help='pheno_covar.fixed.tsv: PHENO1 and covariates')
    p.add_argument('--keep', required=True, help='fixed.keep: the unrelated samples')
    p.add_argument('--tests', required=True, help='tests.<class>.tsv')
    p.add_argument('--thresholds', required=True, help='effective_tests.<class>.tsv')
    p.add_argument('--pheno-col', default='PHENO1')
    p.add_argument('--covar-cols', required=True, help='expanded, for the null fit')
    p.add_argument('--covar-name', required=True, help='plink2 --covar-name, as ASSOC_FIXED')
    p.add_argument('--firth-mode', required=True, help='as ASSOC_FIXED')
    p.add_argument('--plink2', required=True)
    p.add_argument('--alpha', type=float, default=0.05)
    p.add_argument('--draws', type=int, default=10000)
    p.add_argument('--batch', type=int, default=500)
    p.add_argument('--seed', type=int, default=20261023)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--out', required=True)
    p.add_argument('--out-minp', required=True)
    p.add_argument('--out-observed', required=True)
    return p.parse_args()


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f'ERROR: {cmd[0]} failed:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}')


def null_fit(y, X, b=None, tol=1e-10, maxiter=100):
    """Logistic maximum likelihood by Newton-Raphson; returns (coefficients, fitted probabilities)."""
    b = np.zeros(X.shape[1]) if b is None else b.copy()
    for _ in range(maxiter):
        mu = 1 / (1 + np.exp(-X @ b))
        step = np.linalg.solve(X.T @ (X * (mu * (1 - mu))[:, None]), X.T @ (y - mu))
        b += step
        if np.max(np.abs(step)) < tol:
            return b, 1 / (1 + np.exp(-X @ b))
    sys.exit('ERROR: the null logistic model did not converge')


def score_p(y, X, G, G2, b0):
    """Score-test P of every column of G (mean-imputed) on the null model of y; NaN where the variance is 0."""
    b, mu = null_fit(y, X, b0)
    w = mu * (1 - mu)
    u = G.T @ (y - mu)
    A = (X * w[:, None]).T @ G                                   # p x m
    v = G2.T @ w - np.einsum('im,im->m', A, np.linalg.solve(X.T @ (X * w[:, None]), A))
    with np.errstate(invalid='ignore', divide='ignore'):
        z = np.where(v > 1e-12, u / np.sqrt(np.where(v > 1e-12, v, 1)), np.nan)
    return 2 * stats.norm.sf(np.abs(z)), b


def wilson(k, n, z=1.959964):
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def dosage(a, iids, reps):
    """P dosage (samples in iids order x tests), blanks mean-imputed, from the PLINK set."""
    run([a.plink2, '--bfile', a.bfile, '--keep', a.keep, '--extract', 'tests.txt', '--export', 'A',
         '--threads', str(a.threads), '--out', 'dosage'])
    raw = pd.read_csv('dosage.raw', sep='\t', dtype={'FID': str, 'IID': str}, na_values=['NA']).set_index('IID')
    x = {}
    for c in raw.columns[5:]:
        name, allele = c.rsplit('_', 1)
        v = raw[c].to_numpy(dtype=float)
        if allele not in ('A', 'P'):
            sys.exit(f'ERROR: plink2 counted allele {allele!r} for {name}; expected A or P')
        x[name] = 2.0 - v if allele == 'A' else v
    if set(x) != set(reps) or set(raw.index) != set(iids):
        sys.exit('ERROR: exported markers or samples differ from the tests and fixed-model samples')
    G = pd.DataFrame(x, index=raw.index).loc[iids, reps].to_numpy()
    return np.where(np.isnan(G), np.nanmean(G, axis=0), G)


def run_batch(a, tag, pheno, cols):
    """plink2 on one batch of phenotype columns; returns {column: Series of P by marker}."""
    pheno.to_csv(f'{tag}.pheno', sep='\t', index=False)
    run([a.plink2, '--bfile', a.bfile, '--keep', a.keep, '--extract', 'tests.txt',
         '--pheno', f'{tag}.pheno', '--1', '--covar', a.pheno_fixed, '--covar-name', a.covar_name,
         '--glm', 'omit-ref', a.firth_mode, 'hide-covar', '--threads', str(a.threads), '--out', tag])
    out = {}
    for c in cols:
        f = Path(f'{tag}.{c}.glm.logistic')
        if not f.exists():
            sys.exit(f'ERROR: plink2 wrote no {f}')
        out[c] = pd.read_csv(f, sep='\t', usecols=['ID', 'P'], dtype={'ID': str}, na_values=['NA']).set_index('ID').P
        f.unlink()
    Path(f'{tag}.pheno').unlink()
    return out


def summarize(minp, alpha, th, rng):
    thr = float(np.quantile(minp, alpha))
    boot = np.quantile(rng.choice(minp, (1000, len(minp))), alpha, axis=1)
    row = {'p_maxt': thr, 'p_maxt_lo': float(np.percentile(boot, 2.5)), 'p_maxt_hi': float(np.percentile(boot, 97.5)),
           'me_maxt': alpha / thr}
    for t in ('gec', 'bonf'):
        n_hit = int((minp < th[f'p_{t}']).sum())
        lo, hi = wilson(n_hit, len(minp))
        row.update({f'fwer_{t}': n_hit / len(minp), f'fwer_{t}_lo': lo, f'fwer_{t}_hi': hi})
    return row


def main():
    a = parse_args()
    fx = pd.read_csv(a.pheno_fixed, sep='\t', dtype={'#FID': str, 'IID': str})
    keep = pd.read_csv(a.keep, sep=r'\s+', header=None, dtype=str)
    if set(keep[1]) != set(fx.IID) or len(keep) != len(fx):
        sys.exit(f'ERROR: fixed.keep ({len(keep)}) and pheno_covar.fixed.tsv ({len(fx)}) differ')
    cols = a.covar_cols.split(',')
    if fx[[a.pheno_col] + cols].isna().any().any():
        sys.exit('ERROR: missing phenotype or covariate among the fixed-model samples')
    y = fx[a.pheno_col].astype(float).to_numpy()
    if not set(np.unique(y)) <= {0.0, 1.0}:
        sys.exit(f'ERROR: {a.pheno_col} is not 0/1')
    X = np.column_stack([np.ones(len(fx)), fx[cols].astype(float).to_numpy()])
    b_obs, mu = null_fit(y, X)

    tt = pd.read_csv(a.tests, sep='\t', dtype=str, keep_default_na=False)
    reps = tt.loc[tt.test_id == tt.marker, 'marker'].tolist()
    Path('tests.txt').write_text('\n'.join(reps) + '\n')
    th = pd.read_csv(a.thresholds, sep='\t').iloc[0]
    if int(th.tests) != len(reps):
        sys.exit(f'ERROR: {len(reps)} tests in {a.tests}, {int(th.tests)} in {a.thresholds}')
    G = dosage(a, fx.IID.tolist(), reps)
    G2 = G * G
    obs_score, _ = score_p(y, X, G, G2, b_obs)

    rng = np.random.default_rng(a.seed)
    minp = {s: [] for s in STATS}
    ncase, nok, observed = [], {s: [] for s in STATS}, None
    for b0 in range(0, a.draws, a.batch):
        k = min(a.batch, a.draws - b0)
        Y = (rng.random((len(fx), k)) < mu[:, None]).astype(np.int8)
        names = [f'R{b0 + j + 1:05d}' for j in range(k)]
        ph = pd.DataFrame(Y, columns=names)
        if b0 == 0:
            ph['OBSERVED'] = y.astype(np.int8)
        ph.insert(0, 'IID', fx.IID.values)
        ph.insert(0, '#FID', fx['#FID'].values)
        res = run_batch(a, f'b{b0 // a.batch:03d}', ph, list(ph.columns[2:]))
        if b0 == 0:
            observed = res.pop('OBSERVED').reindex(reps)
        for j, c in enumerate(names):
            pw = res[c].reindex(reps).to_numpy()
            ps, _ = score_p(Y[:, j].astype(float), X, G, G2, b_obs)
            for s, pv in (('wald', pw), ('score', ps)):
                minp[s].append(float(np.nanmin(pv)) if np.isfinite(pv).any() else np.nan)
                nok[s].append(int(np.isfinite(pv).sum()))
            ncase.append(int(Y[:, j].sum()))
        print(f'[maxt_boot] {a.case_set} {a.cohort} {a.cls}: {b0 + k} / {a.draws} replicates', flush=True)

    minp = {s: np.array(v) for s, v in minp.items()}
    nok = {s: np.array(v) for s, v in nok.items()}
    ncase = np.array(ncase)
    if any(np.isnan(v).any() for v in minp.values()):
        sys.exit('ERROR: a replicate returned no P for any test')
    row = {'case_set': a.case_set, 'cohort': a.cohort, 'class': a.cls, 'alpha': a.alpha, 'sample_set': 'fixed',
           'samples': len(fx), 'cases': int(y.sum()), 'controls': int(len(y) - y.sum()), 'tests': len(reps),
           'draws': a.draws, 'seed': a.seed, 'null_model': f'logistic {a.pheno_col} ~ 1 + ' + ' + '.join(cols),
           'cases_mean': float(ncase.mean()), 'cases_sd': float(ncase.std(ddof=1)),
           'p_gec': float(th.p_gec), 'p_bonf': float(th.p_bonf)}
    for s in STATS:
        for k, v in summarize(minp[s], a.alpha, th, rng).items():   # p_maxt_lo -> p_maxt_wald_lo
            stem, end = (k[:-3], k[-3:]) if k.endswith(('_lo', '_hi')) else (k, '')
            row[f'{stem}_{s}{end}'] = v
        row.update({f'replicates_with_na_{s}': int((nok[s] < len(reps)).sum()), f'tests_ok_min_{s}': int(nok[s].min())})
    row.update({'p_maxt_mvn': float(th.p_maxt), 'fwer_gec_mvn': float(th.fwer_gec), 'fwer_bonf_mvn': float(th.fwer_bonf),
                'observed_min_p_wald': float(observed.min()), 'observed_min_p_score': float(np.nanmin(obs_score))})
    pd.DataFrame([row]).to_csv(a.out, sep='\t', index=False, float_format='%.6g')
    pd.DataFrame({'replicate': np.arange(1, a.draws + 1), 'cases': ncase,
                  'tests_ok_wald': nok['wald'], 'min_p_wald': minp['wald'],
                  'tests_ok_score': nok['score'], 'min_p_score': minp['score']}) \
        .to_csv(a.out_minp, sep='\t', index=False, float_format='%.6g', compression='gzip')
    pd.DataFrame({'marker': reps, 'p_wald': observed.to_numpy(), 'p_score': obs_score}) \
        .to_csv(a.out_observed, sep='\t', index=False, float_format='%.6g')
    print(f"[maxt_boot] {a.case_set} {a.cohort} {a.cls}: {len(fx)} samples ({int(y.sum())} cases), {len(reps)} tests, "
          f"{a.draws} replicates. Exact threshold: Wald {row['p_maxt_wald']:.3g}, score {row['p_maxt_score']:.3g} "
          f"(multivariate normal {th.p_maxt:.3g}). FWER at GEC: Wald {row['fwer_gec_wald']:.3f}, score "
          f"{row['fwer_gec_score']:.3f} (normal {th.fwer_gec:.3f}); at Bonferroni: Wald {row['fwer_bonf_wald']:.3f}, "
          f"score {row['fwer_bonf_score']:.3f} (normal {th.fwer_bonf:.3f})")


if __name__ == '__main__':
    main()
