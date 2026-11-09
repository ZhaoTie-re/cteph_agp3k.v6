#!/usr/bin/env python3
"""Control-source sensitivity across cohorts (reads results/; run in the report stage).

Every case against the controls of one source (ctrlonly_<source>) and against all sources but one
(ctrlwithout_<source>), each signal with its main model: SAIGE (GWAS signal variants, HLA alleles and
residues; the subset's own null model), rvtest CMC and SKAT-O (Tier 1-2 genes, both variant sets), the fixed
omnibus (residue positions). The five ctrlonly estimates share every case, so they are correlated:
rho_kl = sqrt(m_k m_l / ((m_k + n)(m_l + n))), m = the source's controls, n = the cases (Lin & Sullivan 2009
with the roles of cases and controls swapped). Q, I2, REML and Hartung-Knapp as in platform_stats.py.

  controls.counts.tsv         per cohort x subset: cases, controls (random, fixed, rare-variant sets)
  controls.strata.tsv         per cohort x marker x subset: beta, SE, OR (95% CI), P, allele frequencies
  controls.heterogeneity.tsv  per cohort x marker: pooled, GLS across the five sources, Q, P, I2, RE
  controls.loo.tsv            per cohort x marker x left-out source: beta, SE, P, z observed and expected
                              (beta_pooled / SE without the source), ctrlonly - ctrlwithout and its P
  controls.null_lambda.tsv    per cohort x subset: lambda and share |z| > 3 over the null SNPs
  controls.genes.tsv          per cohort x gene x set x method x subset: P, OR, carriers (cases, controls)
  controls.positions.tsv      per cohort x position: omnibus P per subset
  controls.signals.tsv        the GWAS leads and DQB1:57:V: pooled, the five sources, Q, the weakest LOO
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402
import platform_stats as PS  # noqa: E402
import report_lib as RL  # noqa: E402

SOURCES, ONLY, WITHOUT = RL.SOURCES, RL.SOURCE_ONLY, RL.SOURCE_WITHOUT
SUBSETS = ['pooled'] + ONLY + WITHOUT


def het_rows(sa, markers, n_case, m, c):
    het, loo = [], []
    for vid, mtype in markers:
        if vid not in sa['pooled'].index:
            continue
        p = sa['pooled'].loc[vid]
        row = {'cohort': c, 'variant_id': vid, 'marker_type': mtype, 'beta_pooled': p.beta, 'se_pooled': p.se,
               'p_pooled': p.p}
        b, s, used = [], [], []
        for src, sub in zip(SOURCES, ONLY):
            x = sa[sub].loc[vid] if vid in sa[sub].index else None
            if x is not None and np.isfinite(x.beta) and np.isfinite(x.se) and x.se > 0:
                b.append(x.beta)
                s.append(x.se)
                used.append(src)
        row['sources_used'] = len(used)
        if len(used) >= 2:
            b, s = np.array(b), np.array(s)
            R = np.array([[1.0 if i == j else PS.rho_shared(m[used[i]], m[used[j]], n_case) for j in range(len(used))]
                          for i in range(len(used))])
            mu, se_fe, q = PS.gls(b, R * np.outer(s, s))
            t2, mu_re, se_hk, lo, hi = PS.reml(b, R * np.outer(s, s))
            i2 = max(0.0, (q - (len(used) - 1)) / q) if q > 0 else 0.0
            i2lo, i2hi = PS.i2_ci(q, len(used))
            row.update({'beta_fe': mu, 'se_fe': se_fe, 'Q': q, 'df': len(used) - 1, 'p_Q': stats.chi2.sf(q, len(used) - 1),
                        'I2': i2, 'I2_l95': i2lo, 'I2_u95': i2hi, 'tau2': t2, 'beta_re': mu_re, 'se_re_hk': se_hk,
                        're_l95': lo, 're_u95': hi})
        het.append(row)
        m_all = sum(m.values())
        for src, only, wo in zip(SOURCES, ONLY, WITHOUT):
            x = sa[wo].loc[vid] if vid in sa[wo].index else None
            if x is None or not np.isfinite(x.beta):
                continue
            rec = {'cohort': c, 'variant_id': vid, 'marker_type': mtype, 'left_out': src, 'beta': x.beta, 'se': x.se,
                   'p': x.p, 'n_ctrl': int(x.n_ctrl), 'z_obs': x.beta / x.se, 'z_expected': p.beta / x.se,
                   'beta_pooled': p.beta}
            k_ = sa[only].loc[vid] if vid in sa[only].index else None
            if k_ is not None and np.isfinite(k_.beta):
                r = PS.rho_shared(m[src], m_all - m[src], n_case)
                v = k_.se ** 2 + x.se ** 2 - 2 * r * k_.se * x.se
                d = k_.beta - x.beta
                rec.update({'beta_source_only': k_.beta, 'diff': d, 'se_diff': np.sqrt(v),
                            'p_diff': 2 * stats.norm.sf(abs(d) / np.sqrt(v))})
            loo.append(rec)
    return het, loo


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results', required=True)
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    res, out = Path(a.results), Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    sig = RL.read_signals(res)
    gw = sig[sig.kind == 'gwas']
    nul = pd.read_csv(res / '02.targets' / 'null_snps.tsv', sep='\t', dtype={'variant_id': str})

    counts, strata, het, loo, lam, signals = [], [], [], [], [], []
    for c in K.COHORTS:
        cf = RL.subset_counts_file(res, c, ONLY[0])
        if not cf.exists():
            continue
        cnt = pd.read_csv(cf, sep='\t')
        pooled_cnt = pd.read_csv(RL.subset_counts_file(res, c, 'pooled'), sep='\t')
        cnt = pd.concat([pooled_cnt[pooled_cnt.subset == 'pooled'].assign(source='all'), cnt], ignore_index=True)
        counts.append(cnt)
        sa = RL.read_saige_subsets(res, c, SUBSETS)
        missing = [s for s in SUBSETS if s not in sa]
        if missing:
            sys.exit(f'ABORT: {c}: no SAIGE results for {missing}')
        n_case = int(cnt.set_index('subset').loc['pooled', 'cases_random'])
        m = {s: int(cnt.set_index('subset').loc[f'ctrlonly_{s}', 'controls_random']) for s in SOURCES}

        hla = sorted(i for i in sa['pooled'].index if ('*' in i) or (':' in i and not i.startswith('chr')))
        markers = [(v, 'gwas') for v in sorted(set(gw.variant_id))] + \
                  [(v, 'hla_allele' if '*' in v else 'hla_residue') for v in hla]
        keep = {v for v, _ in markers}
        for sub, d in sa.items():
            x = d[d.index.isin(keep)].reset_index().drop(columns=['subset', 'cohort'], errors='ignore')
            x.insert(0, 'subset', sub)
            x.insert(0, 'cohort', c)
            strata.append(x)
            z = d[d.index.isin(nul.variant_id)].z.dropna()
            lam.append({'cohort': c, 'subset': sub, 'null_snps': len(z),
                        'lambda': float(np.median(z ** 2) / stats.chi2.ppf(0.5, 1)),
                        'share_abs_z_gt3': float((z.abs() > 3).mean())})
        h, lo = het_rows(sa, markers, n_case, m, c)
        het += h
        loo += lo

        hh = pd.DataFrame(h).set_index('variant_id')
        ll = pd.DataFrame(lo)
        for unit, vid in RL.saige_signals(sig):
            if vid not in hh.index:
                continue
            r = hh.loc[vid]
            row = {'cohort': c, 'signal': unit, 'marker': vid, 'or_pooled': np.exp(r.beta_pooled), 'p_pooled': r.p_pooled,
                   'p_Q': r.get('p_Q'), 'I2': r.get('I2'), 'or_re': np.exp(r.get('beta_re', np.nan))}
            for s, only in zip(SOURCES, ONLY):
                if vid in sa[only].index:
                    row[f'or_{s}'] = np.exp(sa[only].loc[vid, 'beta'])
            x = ll[ll.variant_id == vid]
            if len(x):
                ratio = (x.z_obs / x.z_expected)
                w = x.loc[ratio.idxmin()]
                row.update({'weakest_loo_source': w.left_out, 'weakest_loo_z_obs': w.z_obs,
                            'weakest_loo_z_expected': w.z_expected, 'weakest_loo_p': w.p,
                            'min_p_diff': x.p_diff.min() if 'p_diff' in x else np.nan})
            signals.append(row)

    if not counts:
        sys.exit('ABORT: no control-source results found')
    pd.concat(counts).to_csv(out / 'controls.counts.tsv', sep='\t', index=False)
    st = pd.concat(strata, ignore_index=True)
    st['or'] = np.exp(st.beta)
    st['or_l95'] = np.exp(st.beta - 1.959964 * st.se)
    st['or_u95'] = np.exp(st.beta + 1.959964 * st.se)
    st.to_csv(out / 'controls.strata.tsv', sep='\t', index=False)
    pd.DataFrame(het).to_csv(out / 'controls.heterogeneity.tsv', sep='\t', index=False)
    pd.DataFrame(loo).to_csv(out / 'controls.loo.tsv', sep='\t', index=False)
    pd.DataFrame(lam).to_csv(out / 'controls.null_lambda.tsv', sep='\t', index=False)
    pd.DataFrame(signals).to_csv(out / 'controls.signals.tsv', sep='\t', index=False)

    g = []
    for c in K.COHORTS:
        for sub in SUBSETS:
            for f in sorted((res / '09.platform_assoc' / c).glob(f'gene.{sub}.*.tsv')):
                x = pd.read_csv(f, sep='\t')
                x['cohort'], x['subset'] = c, sub
                g.append(x)
    if g:
        keep = ['cohort', 'subset', 'gene', 'stratum', 'method', 'n_var', 'n_informative', 'pvalue', 'beta', 'se', 'or',
                'or_l95', 'or_u95', 'carriers_case', 'carriers_control', 'multi_site_case', 'multi_site_control']
        g = pd.concat(g, ignore_index=True)
        g[[k for k in keep if k in g.columns]].to_csv(out / 'controls.genes.tsv', sep='\t', index=False)

    om, w = RL.omnibus_wide([(sub, RL.omnibus_file(res, sub)) for sub in SUBSETS if RL.omnibus_file(res, sub).exists()])
    if om is not None:
        w.to_csv(out / 'controls.positions.tsv', sep='\t', index=False)
    print(f'[control_report] {len(signals)} headline rows, {len(het)} markers x cohorts -> {out}')


if __name__ == '__main__':
    main()
