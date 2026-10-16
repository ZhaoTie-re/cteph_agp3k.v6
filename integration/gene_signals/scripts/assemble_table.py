#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The long tables behind the workbook, and every cross-check against
#           the upstream tables.
#             genes.long     gene x variant set x cohort: test result (SKAT-O rho
#                            and P, CMC OR (95% CI) and P, BH q, calls), burden
#                            (variants, carriers, cumulative MAC / MAF, multi-site
#                            carriers), the consistency tier of the variant set
#             variants.long  gene x variant x cohort, in the fields of the common-
#                            variant table: rsID, Ensembl 86 gene, snpEff, Open
#                            Targets summary, effect-allele (ALT) frequencies in this
#                            study and the reference panels; per cohort genotypes,
#                            EAF, missingness, HWE P, and whether the cohort's test
#                            used the variant and why not
#             carriers       carriers by case platform and control source
#             scans          per cohort x variant set x statistic: genes tested,
#                            thresholds, calls, lambda_GC
#             lovo           gene x set x cohort x statistic x left-out variant: P
#                            without the variant, delta log10 P, top contributor,
#                            primary variant set (baseline = the published P)
# Inputs  : the step outputs and the upstream tables (see gene_signals.nf)
# Outputs : gene_signals.{genes,variants}.long.tsv, gene_signals.carriers.tsv,
#           gene_signals.scans.tsv, gene_signals.lovo.tsv, checks.json
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gs_common import (csv_list, family_calls, load_json, load_variants, log, num, read_table, read_tsv,  # noqa
                       tiers_from_calls, write_tsv)

METHODS = ('skato', 'cmc')


def primary_sets(G):
    """(set_name, stratum) of each gene's primary variant set(s): the set(s) with a consistency tier; for a gene
    with none, the set(s) with a Bonferroni-significant result; for a gene selected by the sensitivity analysis
    alone, the set(s) with a sensitivity tier."""
    out = set()
    for sn, g in G.groupby('set_name'):
        t = g[g.tier.astype(int) > 0].stratum.unique()
        if not len(t):
            t = g[(g.skato_bonferroni == True) | (g.cmc_bonferroni == True)].stratum.unique()  # noqa: E712 (NaN: untested)
        if not len(t):
            t = g[g.tier_sens.fillna(0).astype(int) > 0].stratum.unique()
        assert len(t), (sn, 'no tier, Bonferroni result or sensitivity tier: why was it selected?')
        out |= {(sn, s) for s in t}
    return out


def lovo_table(paths, G, VL, strata, cohorts):
    """The re-tests, checked against the gene and variant tables: per gene x set x cohort tested, the baseline P is
    the published one and the variants left out are exactly the variants the cohort's test used. Adds
    p_all (the baseline), delta_log10 = log10(P without the variant / P with all), its rank and the top
    contributor (the largest positive delta) per gene x set x cohort x statistic, and the primary-set flag."""
    L = pd.concat([read_tsv(f) for f in paths], ignore_index=True)
    for c in ('k', 'n_var', 'n_left', 'n_var_engine', 'n_poly_engine'):
        L[c] = L[c].astype(int)
    L['p'] = num(L.pvalue)
    gt = G[G.tested].set_index(['set_name', 'stratum', 'cohort'])
    keys = set(map(tuple, L[['set_name', 'stratum', 'cohort']].drop_duplicates().values))
    assert keys == set(gt.index), sorted(keys ^ set(gt.index))[:5]
    for (sn, s, c), g in L.groupby(['set_name', 'stratum', 'cohort']):
        x = gt.loc[(sn, s, c)]
        used = set(VL[(VL.set_name == sn) & (VL.cohort == c) & VL[f'in_set_{s}']].variant_id)
        for m in METHODS:
            gm = g[g.method == m]
            assert set(gm[gm.k > 0].dropped) == used and len(gm[gm.k > 0]) == len(used) == int(x.n_var), (sn, s, c, m)
            b = gm[gm.k == 0]
            assert len(b) == 1 and float(b.p.iloc[0]) == float(x[f'{m}_p']), (sn, s, c, m, b.p.tolist(), x[f'{m}_p'])
    base = L[L.k == 0].set_index(['set_name', 'stratum', 'cohort', 'method']).p
    L['p_all'] = [base.loc[k] for k in zip(L.set_name, L.stratum, L.cohort, L.method)]
    L['delta_log10'] = np.where(L.k > 0, np.log10(L.p) - np.log10(L.p_all), np.nan)
    grp = L[L.k > 0].groupby(['set_name', 'stratum', 'cohort', 'method']).delta_log10
    L['rank'] = grp.rank(ascending=False, method='min')
    L['top'] = (L.k > 0) & (L.delta_log10 > 0) & (L.delta_log10 == grp.transform('max').reindex(L.index))
    prim = primary_sets(G)
    L['primary'] = [(sn, s) in prim for sn, s in zip(L.set_name, L.stratum)]
    order = G.drop_duplicates('set_name').set_index('set_name').gene_order
    L['gene_order'] = L.set_name.map(order).astype(int)
    L['s_rank'] = L.stratum.map({s: i for i, s in enumerate(strata)})
    L['c_rank'] = L.cohort.map({c: i for i, c in enumerate(cohorts)})
    L = L.sort_values(['gene_order', 's_rank', 'c_rank', 'method', 'k']).drop(columns=['s_rank', 'c_rank', 'pvalue'])
    lead = ['gene_order', 'set_name', 'gene_symbol', 'stratum', 'primary', 'cohort', 'method', 'k', 'dropped',
            'n_var', 'n_left', 'n_var_engine', 'n_poly_engine', 'p', 'p_all', 'delta_log10', 'rank', 'top']
    return L[lead + [c for c in L.columns if c not in lead]].reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--genes', required=True)
    ap.add_argument('--variants', required=True)
    ap.add_argument('--selection-json', required=True)
    ap.add_argument('--geno', nargs='+', required=True, help='geno.<cohort>.tsv')
    ap.add_argument('--carriers', nargs='+', required=True, help='carriers.<cohort>.tsv')
    ap.add_argument('--samples', nargs='+', required=True, help='samples.<cohort>.tsv')
    ap.add_argument('--multisite', nargs='+', required=True, help='multisite.<cohort>.tsv (sample IDs)')
    ap.add_argument('--multi-upstream', required=True, help='upstream multi_carriers_all.tsv (cross-check)')
    ap.add_argument('--sample-info', required=True, help='sample sheet (csv): ID_JHRPv6, WGS_Platform, Project, depth')
    ap.add_argument('--v6-info', required=True, help="the project's sample table (xlsx): sex, ages, DVT, mPAP")
    ap.add_argument('--agp3k-agesex', required=True, help='AGP3K age / sex release (tsv): control ages')
    ap.add_argument('--function', required=True)
    ap.add_argument('--tommo', required=True)
    ap.add_argument('--kg', required=True)
    ap.add_argument('--annovar', required=True)
    ap.add_argument('--opentargets', required=True, help='opentargets.tsv (one row per variant)')
    ap.add_argument('--gene-scan', required=True)
    ap.add_argument('--gene-index', required=True)
    ap.add_argument('--scan-qc', required=True)
    ap.add_argument('--denominators', required=True)
    ap.add_argument('--cohort-samples', required=True)
    ap.add_argument('--tiers', required=True)
    ap.add_argument('--rv-qc', nargs='+', required=True, help='rv_qc.<cohort>.json')
    ap.add_argument('--variant-sets', required=True, help='variant_sets.json')
    ap.add_argument('--upstream-manifest', required=True)
    ap.add_argument('--cohorts', required=True)
    ap.add_argument('--cohort-codes', required=True)
    ap.add_argument('--strata', required=True)
    ap.add_argument('--af-cohort', required=True, help='the cohort whose case / control EAF sit beside the references')
    ap.add_argument('--sens-min-num-var', type=int, required=True, help='minimum variants per gene, sensitivity family')
    ap.add_argument('--out-genes', required=True)
    ap.add_argument('--out-variants', required=True)
    ap.add_argument('--out-carriers', required=True)
    ap.add_argument('--out-scans', required=True)
    ap.add_argument('--out-multisite', required=True, help='multi-site carriers by sample ID, with sample information')
    ap.add_argument('--lovo', nargs='+', required=True, help='lovo.<cohort>.<set>.tsv (leave-one-variant-out)')
    ap.add_argument('--out-lovo', required=True)
    ap.add_argument('--out-checks', required=True)
    a = ap.parse_args()
    cohorts, strata = csv_list(a.cohorts), csv_list(a.strata)
    code = dict(zip(cohorts, csv_list(a.cohort_codes)))
    up = load_json(a.upstream_manifest)
    minac = {c: int(up['minac'][c]) for c in cohorts}
    min_num_var = int(up['min_num_var'])
    alpha = float(up['alpha'])
    checks = {}

    genes = read_table(a.genes)
    v = load_variants(a.variants)
    sel = load_json(a.selection_json)

    # ── gene x set x cohort ───────────────────────────────────────────────
    scan = read_tsv(a.gene_scan)
    scan = scan[scan.set_name.isin(genes.set_name) & scan.cohort.isin(cohorts) & scan.stratum.isin(strata)]
    idx = read_tsv(a.gene_index)
    idx = idx[idx.set_name.isin(genes.set_name) & idx.cohort.isin(cohorts) & idx.stratum.isin(strata)]
    den = read_tsv(a.denominators).set_index(['cohort', 'stratum'])
    smp = read_tsv(a.cohort_samples).set_index('cohort')
    tiers = read_tsv(a.tiers)
    car = pd.concat([read_table(f) for f in a.carriers], ignore_index=True)
    tot = car[car.by == 'all'].pivot_table(index=['set_name', 'stratum', 'cohort'], columns='group',
                                           values=['carriers', 'mac', 'multi_site', 'n_var_tested'], aggfunc='first')
    rows = []
    for _, gr in genes.iterrows():
        for s in strata:
            t = tiers[(tiers.set_name == gr.set_name) & (tiers.stratum == s)]
            for c in cohorts:
                rec = {'gene_order': gr.gene_order, 'set_name': gr.set_name, 'gene_symbol': gr.gene_symbol,
                       'chrom': gr.chrom, 'class': gr['class'], 'stratum': s, 'cohort': c,
                       'tier': gr[f'tier_{s}'], 'cohorts_called': t.cohorts_called.iloc[0] if len(t) else '',
                       'n_case': int(smp.loc[c, 'n_case']), 'n_ctrl': int(smp.loc[c, 'n_control']),
                       'bonferroni_threshold': float(den.loc[(c, s), 'bonferroni_threshold']),
                       'n_genes_tested': int(den.loc[(c, s), 'n_genes_mapped'])}
                ix = idx[(idx.set_name == gr.set_name) & (idx.stratum == s) & (idx.cohort == c)]
                rec['tested'] = len(ix) == 1
                if not rec['tested']:
                    rec['not_tested_reason'] = f'fewer than {min_num_var} qualifying variants'
                    rows.append(rec)
                    continue
                rec['not_tested_reason'] = ''
                rec['n_var'] = int(ix.n_var_map.iloc[0])
                k = (gr.set_name, s, c)
                rec.update({'carriers_case': int(tot.loc[k, ('carriers', 'case')]),
                            'carriers_ctrl': int(tot.loc[k, ('carriers', 'ctrl')]),
                            'mac_case': int(tot.loc[k, ('mac', 'case')]), 'mac_ctrl': int(tot.loc[k, ('mac', 'ctrl')]),
                            'multi_case': int(tot.loc[k, ('multi_site', 'case')]),
                            'multi_ctrl': int(tot.loc[k, ('multi_site', 'ctrl')])})
                assert int(tot.loc[k, ('n_var_tested', 'case')]) == rec['n_var'], k
                rec['maf_case'] = rec['mac_case'] / (2 * rec['n_case'])
                rec['maf_ctrl'] = rec['mac_ctrl'] / (2 * rec['n_ctrl'])
                for m in METHODS:
                    x = scan[(scan.set_name == gr.set_name) & (scan.stratum == s) & (scan.cohort == c)
                             & (scan.method == m)]
                    assert len(x) == 1, (k, m, len(x))
                    x = x.iloc[0]
                    # the burden the test saw is the burden counted here
                    assert (int(x.n_carrier_case), int(x.n_carrier_control), int(x.mac_case), int(x.mac_control),
                            int(x.n_var_map)) == (rec['carriers_case'], rec['carriers_ctrl'], rec['mac_case'],
                                                  rec['mac_ctrl'], rec['n_var']), (k, m)
                    rec[f'{m}_p'] = float(x.pvalue)
                    rec[f'{m}_q'] = float(x.fdr_bh)
                    rec[f'{m}_bh'] = x.significant_bh == '1'
                    rec[f'{m}_bonferroni'] = x.significant_bonferroni == '1'
                    assert (float(x.pvalue) < rec['bonferroni_threshold']) == rec[f'{m}_bonferroni'], (k, m)
                    assert (float(x.fdr_bh) < alpha) == rec[f'{m}_bh'], (k, m)
                    if m == 'skato':
                        rec['skato_rho'] = num(pd.Series([x.rho])).iloc[0]
                    else:
                        rec.update({'cmc_or': num(pd.Series([x['or']])).iloc[0],
                                    'cmc_l95': num(pd.Series([x.or_l95])).iloc[0],
                                    'cmc_u95': num(pd.Series([x.or_u95])).iloc[0]})
                rows.append(rec)
    G = pd.DataFrame(rows)
    # the tier of a set is what the calls in the three cohorts give
    for (sn, s), g in G[G.tested].groupby(['set_name', 'stratum']):
        called = g[g.skato_bh | g.cmc_bh].cohort.tolist()
        both = g[g.skato_bh & g.cmc_bh].cohort.tolist()
        tier = 1 if (len(called) == len(cohorts) and both) else 2 if len(called) == len(cohorts) \
            else 3 if len(called) >= 2 else 0
        assert tier == int(g.tier.iloc[0]), (sn, s, tier, g.tier.iloc[0])
    checks['genes'] = {'rows': len(G), 'tested': int(G.tested.sum()), 'not_tested': int((~G.tested).sum())}

    # ── the sensitivity family (genes with >= SensMinNumVar variants): BH q, calls, tier ──
    scan_all = read_tsv(a.gene_scan)
    scan_all = scan_all[scan_all.cohort.isin(cohorts) & scan_all.stratum.isin(strata)]
    sens = family_calls(scan_all, a.sens_min_num_var, alpha)
    t_sens = tiers_from_calls(sens, cohorts)
    sx = sens.set_index(['set_name', 'stratum', 'cohort', 'method'])
    for col in ('in_sens', 'skato_q_sens', 'skato_bh_sens', 'skato_bonferroni_sens', 'cmc_q_sens', 'cmc_bh_sens',
                'cmc_bonferroni_sens', 'bonferroni_threshold_sens', 'n_genes_tested_sens', 'tier_sens'):
        G[col] = np.nan
    G['in_sens'] = False
    for i, r in G.iterrows():
        G.at[i, 'tier_sens'] = t_sens.get((r.set_name, r.stratum), np.nan)
        for m in METHODS:
            k = (r.set_name, r.stratum, r.cohort, m)
            if k not in sx.index:
                continue
            x = sx.loc[k]
            assert r.tested and int(r.n_var) >= a.sens_min_num_var and abs(float(x.p) - r[f'{m}_p']) < 1e-12 * max(
                1, r[f'{m}_p']) + 1e-300, k
            G.at[i, 'in_sens'] = True
            G.at[i, f'{m}_q_sens'] = float(x.q)
            G.at[i, f'{m}_bh_sens'] = bool(x.bh)
            G.at[i, f'{m}_bonferroni_sens'] = bool(x.bonf)
            G.at[i, 'bonferroni_threshold_sens'] = float(x.bonf_thr)
            G.at[i, 'n_genes_tested_sens'] = int(x.family_size)
        # a gene is in the sensitivity family exactly when its cohort's test used >= the sensitivity minimum
        assert G.at[i, 'in_sens'] == (bool(r.tested) and int(r.n_var) >= a.sens_min_num_var), (r.set_name, r.stratum,
                                                                                                  r.cohort)
    gsel = read_table(a.genes).set_index('set_name')
    for (sn, st), g in G.groupby(['set_name', 'stratum']):
        t_file = gsel.loc[sn, f'tier_sens_{st}']
        t_here = g.tier_sens.iloc[0]
        assert (pd.isna(t_file) and not g.in_sens.any()) or int(t_file) == int(t_here), (sn, st, t_file, t_here)
    checks['sensitivity'] = {'min_num_var': a.sens_min_num_var,
                             'tiers': {f'{sn}|{st}': t for (sn, st), t in sorted(t_sens.items()) if t
                                       and sn in set(genes.set_name)}}

    # ── gene x variant x cohort ───────────────────────────────────────────
    geno = pd.concat([read_table(f) for f in a.geno], ignore_index=True)
    samples = {}
    for f in a.samples:
        s_ = read_tsv(f)
        c = next(c for c in cohorts if Path(f).name == f'samples.{c}.tsv')
        samples[c] = set(s_.IID)
    # absent from a cohort's call set for no listed reason: the cohort holds no carrier (the call set keeps
    # variants with a minor allele in its own samples); proven from the cohorts where the variant is present
    carriers_any = {}
    for _, r in geno[geno.in_callset.astype(str) == 'True'].iterrows():
        if isinstance(r.carrier_iids, str) and r.carrier_iids:
            carriers_any.setdefault(r.variant_id, set()).update(r.carrier_iids.split(';'))
    reason = []
    for _, r in geno.iterrows():
        if str(r.in_callset) == 'True' or r.callset_reason != 'not_in_callset':
            reason.append(r.callset_reason if str(r.in_callset) != 'True' else '')
            continue
        cs_ = carriers_any.get(r.variant_id, set())
        assert cs_ and not (cs_ & samples[r.cohort]), (r.variant_id, r.cohort, 'absent but carried')
        reason.append('no_carrier')
    geno['callset_reason'] = reason
    for s in strata:                                     # the set reason carries the call-set reason
        geno[f'set_reason_{s}'] = geno[f'set_reason_{s}'].replace('not_in_callset', 'no_carrier')
    fun = read_table(a.function)
    tom = read_table(a.tommo)
    kg = read_table(a.kg)
    ann = read_table(a.annovar)
    ot = read_table(a.opentargets)
    V = v.merge(fun, on=['set_name', 'variant_id'], validate='1:1', suffixes=('', '_fun'))
    assert (V.impact == V.impact_fun).all() and (V.effect == V.consequence).all()
    V = V.drop(columns=['impact_fun', 'consequence'])
    for t, name in ((tom, 'tommo'), (kg, 'kg'), (ann, 'annovar'), (ot, 'opentargets')):
        assert set(t.variant_id) == set(v.variant_id), name
        V = V.merge(t, on='variant_id', how='left', validate='m:1')
    # rsID as in the common-variant table: the ToMMo record's ID; where it lists several accessions (a
    # multi-allelic position), the one dbSNP 151 assigns to this allele; absent from ToMMo, dbSNP 151's
    tom_rs = V.tommo_rsid.fillna('').astype(str).replace('.', '')
    cands = tom_rs.map(lambda x: sorted({t for t in x.replace(',', ';').split(';') if t.startswith('rs')}))
    db = V.dbsnp_rsid.fillna('').astype(str)
    rs, src, unres = [], [], []
    for cs_, d in zip(cands, db):
        if len(cs_) == 1:
            rs.append(cs_[0]); src.append('ToMMo'); unres.append('')
        elif len(cs_) > 1 and d in cs_:
            rs.append(d); src.append('dbSNP allele'); unres.append('')
        elif len(cs_) > 1:
            rs.append(''); src.append(''); unres.append('/'.join(cs_))
        elif d:
            rs.append(d); src.append('dbSNP 151'); unres.append('')
        else:
            rs.append(''); src.append(''); unres.append('')
    V['rsID'], V['rsID_source'], V['rsID_unresolved'] = rs, src, unres
    one = V[pd.Series([len(c_) == 1 for c_ in cands], index=V.index) & (db != '')]
    checks['rsid'] = {k: int((V.rsID_source == k).sum()) for k in ('ToMMo', 'dbSNP allele', 'dbSNP 151')}
    checks['rsid'].update({'unresolved': int((V.rsID_unresolved != '').sum()), 'none': int((V.rsID == '').sum()),
                           'tommo_vs_dbsnp_disagree': int((one.rsID != one.dbsnp_rsid).sum())})
    VL = V.merge(geno, on='variant_id', validate='m:m')
    assert len(VL) == len(V) * len(cohorts)
    # this study's EAF beside the references: the AfCohort's cases (CTEPH) and controls (AGP3K)
    af = VL[VL.cohort == a.af_cohort].drop_duplicates('variant_id').set_index('variant_id')
    VL['af_cteph'] = VL.variant_id.map(af.case_alt_af)
    VL['af_agp3k'] = VL.variant_id.map(af.ctrl_alt_af)
    for s in strata:
        VL[f'in_set_{s}'] = VL[f'in_set_{s}'].astype(str) == 'True'
    VL['in_callset'] = VL.in_callset.astype(str) == 'True'
    # every gene's tested variants in a cohort are the gene index count
    for (sn, s, c), n in G[G.tested].set_index(['set_name', 'stratum', 'cohort']).n_var.items():
        got = int(VL[(VL.set_name == sn) & (VL.cohort == c)][f'in_set_{s}'].sum())
        assert got == n, (sn, s, c, got, n)
        sub = VL[(VL.set_name == sn) & (VL.cohort == c) & VL[f'in_set_{s}']]
        g = G[(G.set_name == sn) & (G.stratum == s) & (G.cohort == c)].iloc[0]
        assert int(sub.case_mac.sum()) == g.mac_case and int(sub.ctrl_mac.sum()) == g.mac_ctrl, (sn, s, c)
    # ── leave-one-variant-out: the re-tests of each gene x set x cohort ─────
    L = lovo_table(a.lovo, G, VL, strata, cohorts)
    checks['lovo'] = {'rows': len(L), 'tests': int((L.k > 0).sum()),
                      'baselines_equal_to_published': int((L.k == 0).sum()),
                      'primary_sets': {sn: sorted(g.stratum.unique()) for sn, g in L[L.primary].groupby('set_name')},
                      'genes_with_two_variants': sorted({f'{r.set_name}|{r.stratum}|{r.cohort}'
                                                         for r in L[L.n_var == 2].itertuples()}),
                      'top': {f'{r.set_name}|{r.stratum}|{r.cohort}|{r.method}': r.dropped
                              for r in L[L.top].itertuples()}}
    checks['variants'] = {'rows': len(VL), 'by_reason': {s: VL[f'set_reason_{s}'].fillna('').value_counts().to_dict()
                                                         for s in strata},
                          'callset_absent': VL[~VL.in_callset].callset_reason.value_counts().to_dict()}
    cov = V.drop_duplicates('variant_id')
    checks['coverage'] = {'n_variants': len(cov), 'tommo': int(cov.tommo_af.notna().sum()),
                          'kg_eas': int(cov.kg_eas_af.notna().sum()), 'kg_jpt': int(cov.kg_jpt_af.notna().sum()),
                          'gnomad': int(cov.gnomad_af.notna().sum()), 'rsid': int((cov.rsID != '').sum()),
                          'ot_in_index': int(cov.ot_in_index.astype(str).isin(['True', 'true', '1', '1.0']).sum()),
                          'ot_clinvar': int((cov.ot_clinvar_n.fillna(0) > 0).sum()),
                          'ot_gwas': int((cov.ot_gwas_n.fillna(0) > 0).sum()),
                          'ot_qtl': int((cov.ot_qtl_n.fillna(0) > 0).sum()),
                          'ot_e2g': int((cov.ot_e2g_n.fillna(0) > 0).sum())}
    # ToMMo agrees with the controls: rare variants, so the check is on the mean, not a correlation
    fl = VL[(VL.cohort == cohorts[-1]) & VL.in_callset].drop_duplicates('variant_id')
    both = fl.dropna(subset=['tommo_af', 'ctrl_alt_af'])
    checks['tommo_vs_controls'] = {'n': len(both), 'mean_tommo': float(both.tommo_af.mean()),
                                   'mean_controls': float(both.ctrl_alt_af.mean()),
                                   'r': float(np.corrcoef(both.tommo_af, both.ctrl_alt_af)[0, 1]) if len(both) > 2
                                   else None}
    for col in ('tommo_af', 'kg_eas_af', 'kg_jpt_af', 'gnomad_af', 'gnomad_eas_af', 'gnomad_nfe_af', 'case_alt_af',
                'ctrl_alt_af', 'case_hwe_p', 'ctrl_hwe_p'):
        x = VL[col].dropna()
        assert ((x >= 0) & (x <= 1)).all(), col

    # ── carriers, scans, qc, all genes ────────────────────────────────────
    car = car.merge(genes[['set_name', 'gene_symbol', 'gene_order']], on='set_name')
    sq = read_tsv(a.scan_qc)
    sq = sq[sq.cohort.isin(cohorts) & sq.stratum.isin(strata)]
    S = sq.merge(read_tsv(a.denominators)[['cohort', 'stratum', 'n_variants', 'n_genes_dropped_lt_minvar']],
                 on=['cohort', 'stratum'])
    S = S.merge(read_tsv(a.cohort_samples), on='cohort')
    assert len(S) == len(cohorts) * len(strata) * len(METHODS), len(S)
    # the genes behind each lambda_GC (upstream: all genes with a P = n_var_map >= MinNumVar, and n_var_map >= 3),
    # and lambda_GC recomputed on them as upstream does (median of the 1-df chi-square of P / its null median)
    from statistics import NormalDist, median
    full = read_tsv(a.gene_scan)
    full = full[full.cohort.isin(cohorts) & full.stratum.isin(strata)]
    full['p'] = num(full.pvalue)
    full['nv'] = num(full.n_var_map)
    lam = []
    for (c, s, m), x in full.groupby(['cohort', 'stratum', 'method']):
        x = x[x.p.notna()]
        rec = {'cohort': c, 'stratum': s, 'method': m, 'n_genes_lambda': len(x),
               'n_genes_lambda_ge3': int((x.nv >= 3).sum()), 'n_genes_lambda_2': int((x.nv == 2).sum())}
        assert (x.nv >= min_num_var).all() and rec['n_genes_lambda'] == rec['n_genes_lambda_2'] + rec['n_genes_lambda_ge3']
        for key, sub in (('lambda_gc', x), ('lambda_gc_nvar_ge3', x[x.nv >= 3]), ('lambda_gc_nvar2', x[x.nv == 2])):
            chi = [NormalDist().inv_cdf(max(p, 1e-300) / 2) ** 2 for p in sub.p]
            rec[f'{key}_recomputed'] = median(chi) / 0.4549364231195729
        lam.append(rec)
    S = S.merge(pd.DataFrame(lam), on=['cohort', 'stratum', 'method'], validate='1:1')
    for col in ('lambda_gc', 'lambda_gc_nvar_ge3', 'lambda_gc_nvar2'):
        assert (abs(num(S[col]) - S[f'{col}_recomputed']) < 5e-4).all(), (col, S[[col, f'{col}_recomputed']])
    assert (S.n_genes_lambda == num(S.n_genes_returned)).all()
    # the families per scan: the primary recomputed (= upstream), and the sensitivity family
    prim = family_calls(scan_all, min_num_var, alpha)
    fam = []
    for (c, s_, m), x in prim.groupby(['cohort', 'stratum', 'method']):
        y = sens[(sens.cohort == c) & (sens.stratum == s_) & (sens.method == m)]
        fam.append({'cohort': c, 'stratum': s_, 'method': m, 'n_bh_recomputed': int(x.bh.sum()),
                    'n_bonferroni_recomputed': int(x.bonf.sum()),
                    'n_genes_sens': len(y), 'bonferroni_threshold_sens': float(y.bonf_thr.iloc[0]),
                    'n_significant_bh_sens': int(y.bh.sum()), 'n_significant_bonferroni_sens': int(y.bonf.sum()),
                    'bh_threshold_sens': float(y[y.bh].p.max()) if y.bh.any() else np.nan,
                    'min_p_sens': float(y.p.min())})
    S = S.merge(pd.DataFrame(fam), on=['cohort', 'stratum', 'method'], validate='1:1')
    assert (S.n_bh_recomputed == num(S.n_significant_bh)).all() and \
        (S.n_bonferroni_recomputed == num(S.n_significant_bonferroni)).all(), 'primary family differs from upstream'
    # the sensitivity family is the >= 3-variant lambda_GC family
    assert (S.n_genes_sens == S.n_genes_lambda_ge3).all() or a.sens_min_num_var != 3
    # the rare-variant QC behind the tests: its analysed samples and variants are the tests' own
    vs = load_json(a.variant_sets)
    for f in a.rv_qc:
        q = load_json(f)
        c = q['cohort']
        if c not in cohorts:
            continue
        k = int(q['minac']['recommended_minac'])
        assert k == minac[c], (c, k, minac[c])
        assert (q['n_kept']['case'], q['n_kept']['control']) == (int(smp.loc[c, 'n_case']), int(smp.loc[c, 'n_control'])), c
        assert q['sweep'][str(k)] == vs[c]['analysed'], (c, q['sweep'][str(k)], vs[c]['analysed'])
        for s in strata:
            assert vs[f'{c}|{s}']['variants'] == int(den.loc[(c, s), 'n_variants']), (c, s)
            assert vs[f'{c}|{s}']['genes'] == int(den.loc[(c, s), 'n_genes_mapped']), (c, s)
        checks.setdefault('qc', {})[c] = {'samples': q['n_kept'], 'analysed_variants': q['sweep'][str(k)], 'minac': k,
                                          'gap_protected_rescued': q['variant_qc']['n_protected_rescued']}
    assert set(checks['qc']) == set(cohorts)

    # ── multi-site carriers: equal to the upstream records, then anonymous ────────
    ms = pd.concat([read_tsv(f) for f in a.multisite], ignore_index=True)
    upm = read_tsv(a.multi_upstream)
    upm = upm[upm.cohort.isin(cohorts) & upm.stratum.isin(strata) & upm.set_name.isin(genes.set_name)]
    def keyset(d, ids, gts):
        return {(r.cohort, r.stratum, r.set_name, r.IID if 'IID' in d else r.iid, r.group,
                 tuple(sorted(zip(getattr(r, ids).split(';'), getattr(r, gts).split(';'))))) for r in d.itertuples()}
    upm = upm.assign(group=upm.group.replace({'control': 'ctrl'}))
    mine, theirs = keyset(ms, 'variant_ids', 'dosages'), keyset(upm, 'variant_ids', 'genotypes')
    assert mine == theirs, ('multi-site carriers differ from upstream', len(mine), len(theirs),
                            sorted(mine ^ theirs)[:3])
    # one row per gene x sample: the variants carried, the sets and cohorts counting the sample, and its
    # sequencing and clinical information (the deliverable lists the samples by ID, at the project's request)
    vpos = v.drop_duplicates(['set_name', 'variant_id']).set_index(['set_name', 'variant_id'])
    rows = []
    for (sn, iid), x in ms.groupby(['set_name', 'IID'], sort=False):
        dose = {}
        for _, r in x.iterrows():
            for vid, d in zip(r.variant_ids.split(';'), r.dosages.split(';')):
                assert dose.get(vid, d) == d, (sn, iid, vid, 'dosage differs between cohorts or sets')
                dose[vid] = d
        vids = sorted(dose, key=lambda q: int(vpos.loc[(sn, q), 'pos']))
        rec = {'set_name': sn, 'IID': iid, 'group': x.group.iloc[0], 'n_variants_carried': len(vids),
               'variant_ids': ';'.join(vids), 'dosages': ';'.join(dose[q] for q in vids),
               'impacts': ';'.join(vpos.loc[(sn, q), 'impact'] for q in vids)}
        for st in strata:
            cs = set(x[x.stratum == st].cohort)
            rec[f'cohorts_{st}'] = ''.join(code[c] for c in cohorts if c in cs)
        rows.append(rec)
    MS = pd.DataFrame(rows)
    # sample information: the sample sheet (platform, project, depth), the project's sample table (sex, case ages,
    # DVT, mPAP) and, for control ages, the AGP3K age / sex release
    cram = pd.read_csv(a.sample_info, dtype=str).set_index('ID_JHRPv6')
    assert cram.index.is_unique
    info = pd.read_excel(a.v6_info, dtype=str)
    info = info[info.Flag_JHRPv6.astype(str).str.lower().isin(['true', '1'])].set_index('ID_JHRPv6')
    assert info.index.is_unique
    rel = read_tsv(a.agp3k_agesex).set_index('ID')
    for col, f in (('project', lambda i: cram.loc[i, 'Project']), ('platform', lambda i: cram.loc[i, 'WGS_Platform']),
                   ('target_depth', lambda i: cram.loc[i, 'Target_Depth']),
                   ('observed_depth', lambda i: float(cram.loc[i, 'Observed_Depth'])),
                   ('sex', lambda i: info.loc[i, 'Sex'])):
        MS[col] = [f(i) for i in MS.IID]
    is_case = MS.group == 'case'
    assert (info.loc[MS.IID, 'Outcome'].to_numpy() == np.where(is_case, 'PH', 'AGP3K')).all(), 'case / control'
    age = []
    for i, c_ in zip(MS.IID, is_case):
        src = info.loc[i, 'Age_at_DNA_Collection'] if c_ else rel.loc[i, 'age']
        age.append(num(pd.Series([src])).iloc[0])
    MS['age_at_dna'] = age
    for col, src in (('age_at_diagnosis', 'Age_at_Diagnosis'), ('dvt', 'DVT_Presence'), ('max_mpap', 'Max_mPAP')):
        MS[col] = [info.loc[i, src] if c_ else np.nan for i, c_ in zip(MS.IID, is_case)]
    MS['age_at_diagnosis'] = num(MS.age_at_diagnosis)
    MS['max_mpap'] = num(MS.max_mpap)
    # the sample sheets agree with the genotypes' own record of platform and source
    src_ms = ms.drop_duplicates('IID').set_index('IID')
    assert (src_ms.loc[MS.IID, 'platform'].to_numpy() == MS.platform.to_numpy()).all()
    assert (src_ms.loc[MS.IID, 'source'].to_numpy() == MS.project.to_numpy()).all()
    MS = MS.merge(genes[['set_name', 'gene_symbol', 'gene_order']], on='set_name')
    MS = MS.sort_values(['gene_order', 'group', 'IID']).reset_index(drop=True)          # cases first
    # every counted multi-site carrier is on a row, in each cohort and set that counted it
    for (sn, s_, c), g in G[G.tested].groupby(['set_name', 'stratum', 'cohort']):
        x = MS[(MS.set_name == sn) & MS[f'cohorts_{s_}'].str.contains(code[c])]
        assert (int((x.group == 'case').sum()), int((x.group == 'ctrl').sum())) == \
            (int(g.multi_case.iloc[0]), int(g.multi_ctrl.iloc[0])), (sn, s_, c)
    checks['multisite'] = {'rows': len(MS), 'samples': int(MS.IID.nunique()),
                           'equal_to_upstream_records': len(mine)}

    checks['selection'] = sel
    checks['all_passed'] = True
    write_tsv(G, a.out_genes)
    keep = [c for c in VL.columns if c != 'carrier_iids']               # sample IDs stay out of the deliverables
    write_tsv(VL[keep], a.out_variants)
    write_tsv(car, a.out_carriers)
    write_tsv(S, a.out_scans)
    write_tsv(MS, a.out_multisite)
    write_tsv(L, a.out_lovo)
    json.dump(checks, open(a.out_checks, 'w'), indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x))
    log(f'assembled: {len(G)} gene x set x cohort rows ({checks["genes"]["tested"]} tested), {len(VL)} variant x '
        f'cohort rows, {len(car)} carrier rows, {len(S)} scans; rsID '
        f'{checks["rsid"]}; coverage {checks["coverage"]}')


if __name__ == '__main__':
    main()
