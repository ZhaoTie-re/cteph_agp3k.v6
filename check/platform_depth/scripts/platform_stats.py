#!/usr/bin/env python3
"""Platform-sensitivity statistics of one cohort (strata share all controls; user decision 2026-10-30).

For each marker (GWAS signal variants, HLA alleles and residues, null SNPs) with per-subset SAIGE results
(beta = BETA, se = 1/sqrt(var), consistent with the primary no-SPA P):

(a) heterogeneity across the case platforms k = 1..K (strata: platform-k cases + all controls)
    Sigma_kl = rho_kl se_k se_l, rho_kl = sqrt(n_k n_l / ((n_k + n0)(n_l + n0)))   (Lin & Sullivan 2009)
    GLS: beta_FE = 1'S^-1 b / 1'S^-1 1;  Q = (b - beta_FE 1)' S^-1 (b - beta_FE 1) ~ chi2(K-1)
    I2 = max(0, (Q - df) / Q), CI by the test-based method of Higgins & Thompson (2002);
    P of Q also empirically, against the null SNPs of the same MAF bin and type (SNV / indel);
    power to see a stratum without effect: Phi(|beta_pooled| / se_k - 1.96) + Phi(-|beta_pooled| / se_k - 1.96)
(b) meta-analysis vs pooled: beta_FE (with shared controls its weights follow the case counts, so it
    tracks the pooled estimate: a check), and random effects: tau2 by REML with V = Sigma + tau2 I,
    Hartung-Knapp CI on t(K-1); against the pooled (published) estimate
(c) leave-one-platform-out: beta, se, P without platform k; z expected if the effect were unchanged,
    beta_pooled / se_(-k); the difference beta_k - beta_(-k), its variance with the shared controls
case-only check: within unrelated cases, multinomial logistic platform ~ marker + sex + 10 PCs, LRT
    with K-1 df (equal ORs across platforms <=> equal marker distribution across case platforms when
    the controls are shared; Begg & Zhang 1994). Genes: carriers x platform, permutation chi-square.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bedio as B  # noqa: E402
import common as K  # noqa: E402

STRATA = [f'only_{K.PLATFORM_SLUG[p]}' for p in K.CASE_PLATFORMS]
LOPO = [f'without_{K.PLATFORM_SLUG[p]}' for p in K.CASE_PLATFORMS]
MAF_BINS = [0.01, 0.02, 0.05, 0.10, 0.20, 0.50]


def i2_ci(q, k):
    df = k - 1
    if q <= 0 or df < 1:
        return np.nan, np.nan
    if q > k:
        se = 0.5 * (np.log(q) - np.log(df)) / (np.sqrt(2 * q) - np.sqrt(2 * k - 3))
    else:
        se = np.sqrt(1 / (2 * (k - 2)) * (1 - 1 / (3 * (k - 2) ** 2))) if k > 2 else np.nan
    h = np.sqrt(q / df)
    lo, hi = np.exp(np.log(h) - 1.959964 * se), np.exp(np.log(h) + 1.959964 * se)
    f = lambda x: max(0.0, (x ** 2 - 1) / x ** 2)
    return f(lo), f(hi)


def gls(b, S):
    one = np.ones(len(b))
    Si = np.linalg.inv(S)
    w = one @ Si @ one
    mu = (one @ Si @ b) / w
    r = b - mu
    return mu, np.sqrt(1 / w), float(r @ Si @ r)


def reml(b, S):
    k = len(b)
    one = np.ones(k)

    def nll(t2):
        V = S + t2 * np.eye(k)
        Vi = np.linalg.inv(V)
        w = one @ Vi @ one
        mu = (one @ Vi @ b) / w
        r = b - mu
        return 0.5 * (np.linalg.slogdet(V)[1] + np.log(w) + r @ Vi @ r)

    hi = max(10 * np.var(b), 1e-8)
    t2 = optimize.minimize_scalar(nll, bounds=(0, hi), method='bounded').x
    if nll(0.0) <= nll(t2):
        t2 = 0.0
    V = S + t2 * np.eye(k)
    Vi = np.linalg.inv(V)
    w = one @ Vi @ one
    mu = (one @ Vi @ b) / w
    r = b - mu
    q_hk = (r @ Vi @ r) / (k - 1)
    se_hk = np.sqrt(q_hk / w)
    t = stats.t.ppf(0.975, k - 1)
    return t2, mu, se_hk, mu - t * se_hk, mu + t * se_hk


def rho_shared(n_k, n_l, n0):
    return np.sqrt(n_k * n_l / ((n_k + n0) * (n_l + n0)))


def maf_bin(m):
    for lo, hi in zip(MAF_BINS[:-1], MAF_BINS[1:]):
        if lo <= m < hi or (hi == 0.5 and m == 0.5):
            return f'{lo}-{hi}'
    return 'below'


def load_saige(specs):
    d = {}
    for spec in specs:
        sub, path = spec.split('=', 1)
        x = pd.read_csv(path, sep='\t', dtype={'variant_id': str})
        d[sub] = x.set_index('variant_id')
    return d


def heterogeneity(sa, counts, markers, null_ref=None, null_bins=None):
    n0 = int(counts.loc['pooled', 'controls_random'])
    nk = {s: int(counts.loc[s, 'cases_random']) for s in STRATA}
    het, lopo = [], []
    for vid, mtype in markers:
        pooled = sa['pooled'].loc[vid] if vid in sa['pooled'].index else None
        if pooled is None or not np.isfinite(pooled.beta):
            continue
        b, s, used = [], [], []
        for st in STRATA:
            x = sa[st].loc[vid] if vid in sa[st].index else None
            if x is not None and np.isfinite(x.beta) and np.isfinite(x.se) and x.se > 0:
                b.append(x.beta)
                s.append(x.se)
                used.append(st)
        row = {'variant_id': vid, 'marker_type': mtype, 'beta_pooled': pooled.beta, 'se_pooled': pooled.se,
               'p_pooled': pooled.p, 'af_pooled': (pooled.af_case * pooled.n_case + pooled.af_ctrl * pooled.n_ctrl) /
               (pooled.n_case + pooled.n_ctrl), 'strata_used': len(used)}
        maf = min(row['af_pooled'], 1 - row['af_pooled'])
        row['maf_bin'] = maf_bin(maf)
        if len(used) >= 2:
            b, s = np.array(b), np.array(s)
            R = np.array([[1.0 if i == j else rho_shared(nk[used[i]], nk[used[j]], n0) for j in range(len(used))]
                          for i in range(len(used))])
            S = R * np.outer(s, s)
            mu, se_fe, q = gls(b, S)
            kk = len(used)
            t2, mu_re, se_hk, lo, hi = reml(b, S)
            i2 = max(0.0, (q - (kk - 1)) / q) if q > 0 else 0.0
            i2lo, i2hi = i2_ci(q, kk)
            row.update({'beta_fe': mu, 'se_fe': se_fe, 'Q': q, 'df': kk - 1, 'p_Q': stats.chi2.sf(q, kk - 1),
                        'I2': i2, 'I2_l95': i2lo, 'I2_u95': i2hi, 'tau2': t2, 'beta_re': mu_re, 'se_re_hk': se_hk,
                        're_l95': lo, 're_u95': hi})
            for st, bb, ss in zip(used, b, s):
                z = abs(pooled.beta) / ss
                row[f'power_null_{st}'] = stats.norm.cdf(z - 1.959964) + stats.norm.cdf(-z - 1.959964)
            if null_ref is not None and mtype != 'null':
                ref = null_ref[(null_ref.maf_bin == row['maf_bin']) &
                               (null_ref.is_indel == null_bins.get(vid, False))]
                row['p_Q_empirical'] = (1 + (ref.Q >= q).sum()) / (1 + len(ref)) if len(ref) else np.nan
                row['null_reference_n'] = len(ref)
        het.append(row)
        for st, lo_name in zip(STRATA, LOPO):
            x = sa[lo_name].loc[vid] if vid in sa[lo_name].index else None
            k_ = sa[st].loc[vid] if vid in sa[st].index else None
            if x is None or not np.isfinite(x.beta):
                continue
            n_k = nk[st]
            n_mk = int(counts.loc[lo_name, 'cases_random'])
            rec = {'variant_id': vid, 'marker_type': mtype, 'left_out': st.replace('only_', ''),
                   'beta': x.beta, 'se': x.se, 'p': x.p, 'p_spa': x.p_spa, 'n_case': int(x.n_case),
                   'z_obs': x.beta / x.se, 'z_expected': pooled.beta / x.se, 'beta_pooled': pooled.beta}
            if k_ is not None and np.isfinite(k_.beta):
                r = rho_shared(n_k, n_mk, n0)
                v = k_.se ** 2 + x.se ** 2 - 2 * r * k_.se * x.se
                d = k_.beta - x.beta
                rec.update({'beta_left_out_platform': k_.beta, 'diff': d, 'se_diff': np.sqrt(v),
                            'p_diff': 2 * stats.norm.sf(abs(d) / np.sqrt(v))})
            lopo.append(rec)
    return pd.DataFrame(het), pd.DataFrame(lopo)


def strata_table(sa, markers):
    rows = []
    keep = set(v for v, _ in markers)
    for sub, d in sa.items():
        x = d[d.index.isin(keep)].reset_index().drop(columns='subset', errors='ignore')
        x.insert(0, 'subset', sub)
        rows.append(x)
    t = pd.concat(rows, ignore_index=True)
    t['or'] = np.exp(t.beta)
    t['or_l95'] = np.exp(t.beta - 1.959964 * t.se)
    t['or_u95'] = np.exp(t.beta + 1.959964 * t.se)
    return t


def case_only(geno, cases, platform, covar, perm_markers=(), n_perm=20000, seed=20261030):
    import statsmodels.api as sm
    y = pd.Categorical(platform.loc[cases], categories=K.CASE_PLATFORMS).codes
    X0 = sm.add_constant(covar.loc[cases].astype(float).values, has_constant='add')
    out = []
    rng = np.random.default_rng(seed)
    for name, g in geno.items():
        gv = g.loc[cases].astype(float).values
        ok = ~np.isnan(gv)
        rec = {'marker': name, 'n': int(ok.sum()), 'carriers_or_copies': float(np.nansum(gv))}
        try:
            m0 = sm.MNLogit(y[ok], X0[ok]).fit(disp=0, maxiter=200)
            m1 = sm.MNLogit(y[ok], np.column_stack([X0[ok], gv[ok]])).fit(disp=0, maxiter=200)
            lrt = 2 * (m1.llf - m0.llf)
            rec.update({'lrt': lrt, 'df': len(K.CASE_PLATFORMS) - 1, 'p_lrt': stats.chi2.sf(lrt, len(K.CASE_PLATFORMS) - 1),
                        'converged': bool(m0.mle_retvals.get('converged', True) and m1.mle_retvals.get('converged', True))})
        except Exception as e:      # separation with very few carriers
            rec.update({'lrt': np.nan, 'df': len(K.CASE_PLATFORMS) - 1, 'p_lrt': np.nan, 'converged': False,
                        'fit_error': type(e).__name__})
        if name in perm_markers:
            c = (gv[ok] > 0).astype(int)
            yy = y[ok]
            tab = lambda yv: np.array([[np.sum((yv == j) & (c == i)) for j in range(len(K.CASE_PLATFORMS))] for i in (0, 1)])

            def chi2(t):
                e = t.sum(1, keepdims=True) * t.sum(0, keepdims=True) / t.sum()
                with np.errstate(divide='ignore', invalid='ignore'):
                    return np.nansum(np.where(e > 0, (t - e) ** 2 / e, 0))
            obs = chi2(tab(yy))
            ge = sum(chi2(tab(rng.permutation(yy))) >= obs - 1e-12 for _ in range(n_perm))
            rec.update({'carriers_by_platform': ','.join(str(int(v)) for v in tab(yy)[1]),
                        'chi2_perm_stat': obs, 'p_perm': (ge + 1) / (n_perm + 1), 'n_perm': n_perm})
        out.append(rec)
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohort', required=True)
    ap.add_argument('--counts', required=True)
    ap.add_argument('--saige', nargs='+', required=True, help='<subset>=<collected SAIGE results (all markers)>')
    ap.add_argument('--null', required=True, help='null_snps.tsv')
    ap.add_argument('--signals', required=True, help='signal_variants.tsv')
    ap.add_argument('--samples', required=True)
    ap.add_argument('--pheno-fixed', required=True, help='pooled pheno_covar.fixed.tsv (unrelated)')
    ap.add_argument('--gwas-bfile', required=True)
    ap.add_argument('--hla-matrices', nargs=2, required=True, help='allele and residue dosage of the cohort')
    ap.add_argument('--carriers-bfile', required=True)
    ap.add_argument('--covar-cols', required=True)
    ap.add_argument('--out-prefix', required=True)
    a = ap.parse_args()

    counts = pd.read_csv(a.counts, sep='\t').set_index('subset')
    sa = load_saige(a.saige)
    for sub in ['pooled'] + STRATA + LOPO:
        if sub not in sa:
            sys.exit(f'ABORT: no SAIGE results for subset {sub}')
    nul = pd.read_csv(a.null, sep='\t', dtype={'variant_id': str})
    sig = pd.read_csv(a.signals, sep='\t', dtype=str)
    gw = sorted(set(sig[sig.kind == 'gwas'].variant_id))
    pooled = sa['pooled']
    hla = sorted(i for i in pooled.index if ':' in i and not i.startswith('chr') or '*' in i)
    markers = [(v, 'gwas') for v in gw] + [(v, 'hla_allele' if '*' in v else 'hla_residue') for v in hla]

    # null SNPs first: their Q is the empirical reference
    nmark = [(v, 'null') for v in nul.variant_id if v in pooled.index]
    nh, _ = heterogeneity(sa, counts, nmark)
    nh = nh.merge(nul[['variant_id', 'is_indel']], on='variant_id')
    nh['maf_bin'] = nh.maf_bin.astype(str)
    is_indel = {v: (len(v.split(':')[2]) != 1 or len(v.split(':')[3]) != 1) for v in gw}
    het, lopo = heterogeneity(sa, counts, markers, null_ref=nh, null_bins=is_indel)
    strata_table(sa, markers).to_csv(f'{a.out_prefix}.strata.tsv', sep='\t', index=False)
    het.to_csv(f'{a.out_prefix}.heterogeneity.tsv', sep='\t', index=False)
    lopo.to_csv(f'{a.out_prefix}.lopo.tsv', sep='\t', index=False)

    # the shared-control correlation: analytic vs the null SNPs' z
    n0 = int(counts.loc['pooled', 'controls_random'])
    rows = []
    for i, s1 in enumerate(STRATA):
        for s2 in STRATA[i + 1:]:
            z1 = sa[s1].loc[nh.variant_id, 'z']
            z2 = sa[s2].loc[nh.variant_id, 'z']
            ok = np.isfinite(z1.values) & np.isfinite(z2.values)
            rows.append({'stratum_1': s1, 'stratum_2': s2,
                         'rho_analytic': rho_shared(counts.loc[s1, 'cases_random'], counts.loc[s2, 'cases_random'], n0),
                         'r_null_snps': float(np.corrcoef(z1.values[ok], z2.values[ok])[0, 1]), 'n_snps': int(ok.sum())})
    pd.DataFrame(rows).to_csv(f'{a.out_prefix}.rho.tsv', sep='\t', index=False)

    # null calibration: per stratum lambda and |z| > 3, by MAF bin and type; Q of the null SNPs
    cal = []
    for sub in ['pooled'] + STRATA + LOPO:
        d = sa[sub].loc[sa[sub].index.isin(nul.variant_id)].join(nul.set_index('variant_id')[['maf_bin', 'is_indel']])
        for (b, ind), g in d.groupby(['maf_bin', 'is_indel']):
            z = g.z.dropna()
            cal.append({'subset': sub, 'maf_bin': b, 'is_indel': ind, 'n': len(z),
                        'lambda': float(np.median(z ** 2) / stats.chi2.ppf(0.5, 1)), 'share_abs_z_gt3': float((z.abs() > 3).mean())})
    cal = pd.DataFrame(cal)
    qref = nh.groupby(['maf_bin', 'is_indel']).agg(n=('Q', 'size'), q_median=('Q', 'median'),
                                                   share_p_Q_lt_0_05=('p_Q', lambda p: float((p < 0.05).mean())),
                                                   share_p_Q_lt_0_001=('p_Q', lambda p: float((p < 0.001).mean()))).reset_index()
    cal.to_csv(f'{a.out_prefix}.null_calibration.tsv', sep='\t', index=False)
    qref.to_csv(f'{a.out_prefix}.null_Q.tsv', sep='\t', index=False)

    # case-only check
    st = K.read_tsv(a.samples).set_index('sample_id')
    ph = pd.read_csv(a.pheno_fixed, sep='\t', dtype=str).set_index('IID')
    cases = [s for s in ph.index if ph.loc[s, 'PHENO1'] == '1']
    covar = ph[a.covar_cols.split(',')]
    geno = {}
    fam = B.read_fam(a.gwas_bfile)
    bim = B.read_bim(a.gwas_bfile)
    want = set(gw)
    idx = [i for i, f in enumerate(bim) if f[1] in want]
    codes = B.unpack(B.read_rows(a.gwas_bfile, idx, len(fam)), len(fam))
    for k, i in enumerate(idx):
        f = bim[i]
        _, _, ref, alt = f[1].split(':', 3)
        c = codes[k].astype(float)
        if (f[4], f[5]) == (alt, ref):
            v = np.where(c == B.HOM5, 2, np.where(c == B.HET, 1, np.where(c == B.HOM6, 0, np.nan)))
        else:
            v = np.where(c == B.HOM5, 0, np.where(c == B.HET, 1, np.where(c == B.HOM6, 2, np.nan)))
        geno[f[1]] = pd.Series(v, index=fam)
    for path in a.hla_matrices:
        m = pd.read_csv(path, sep='\t', dtype={'sample_id': str}).set_index('sample_id')
        for c in m.columns:
            if c in set(hla):
                geno[c] = m[c].astype(float)
    cf = B.read_fam(a.carriers_bfile)
    cb = B.read_bim(a.carriers_bfile)
    cc = B.unpack(B.read_rows(a.carriers_bfile, list(range(len(cb))), len(cf)), len(cf))
    gene_markers = []
    for k, f in enumerate(cb):
        geno[f[1]] = pd.Series((cc[k] == B.HET).astype(float), index=cf)
        gene_markers.append(f[1])
    cases_g = [s for s in cases if all(s in g.index for g in geno.values())]
    co = case_only(geno, cases_g, st.platform, covar, perm_markers=set(gene_markers))
    co.to_csv(f'{a.out_prefix}.case_only.tsv', sep='\t', index=False)

    summ = {'cohort': a.cohort, 'markers': len(markers), 'null_snps': len(nh),
            'gwas_signals_Q_p_lt_0_05': int((het[het.marker_type == 'gwas'].p_Q < 0.05).sum()),
            'hla_markers_Q_p_lt_0_05': int((het[het.marker_type != 'gwas'].p_Q < 0.05).sum()),
            'null_share_p_Q_lt_0_05': float((nh.p_Q < 0.05).mean())}
    Path(f'{a.out_prefix}.summary.json').write_text(json.dumps(summ, indent=2) + '\n')
    print(json.dumps(summ, indent=2))


if __name__ == '__main__':
    main()
