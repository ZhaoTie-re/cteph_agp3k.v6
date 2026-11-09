#!/usr/bin/env python3
"""Array evidence per signal and cohort, from the per-cohort re-run of check/array_crosscheck (array_cohorts.py).

The unrelated genomes with SNP arrays are DNBSeq-T7 cases (ASA24v1) and HiSeq X controls (JSA24v1): narrow
mainland 1,398 (279 cases, 1,119 controls), intermediate 1,460 (289, 1,171), full 1,683 (299, 1,384). NovaSeq and
DNBSeq-G400RS genomes have no arrays. Each comparison is WGS against the array in the same people and at the same
marker and allele; array_cohorts.py re-runs array_crosscheck's steps in each cohort (gate G10: full mainland
reproduces array_crosscheck's own tables).

One question (user 2026-11-05): does the signal look the same when the WGS genotypes are replaced by array data in
the same people? Each signal takes ONE marker, in this order, and cases and controls are read from the same kind of
data:
  typed             the lead itself typed on the arrays (both groups)
  imputed           the lead in the imputed arrays (Minimac4 v4.1.4, AGP3K panel) with imputation R2 >= 0.8; the R2
                    is the imputation's own (one value; identical in the three cohorts); DQB1:57:V: DQB1 imputed with
                    HIBAG (array_hibag.py)
  LD proxy          a typed proxy of the lead with r2 >= 0.8 (r2 in the cohort's controls)
  none              not checkable
Odds ratios are per the ALT allele of the lead (DQB1:57:V: per V; an LD proxy: per its A1).
Each marker also carries its genotype concordance with WGS in cases and in controls (user 2026-11-05): typed, the
array genotype (04.compare discordance counts); imputed, the best-guess genotype (06.imputed/concordance.tsv);
DQB1:57:V, the HIBAG best-guess 57 V count against the HLA-HD one.
The association is re-run on that marker in the same people (sex + 10 PCs) and compared with WGS:
  contradicted = opposite direction or |z WGS - z array| > 3; otherwise consistent (user 2026-11-05: array P < 0.05
  is not required; the array P and both z are kept as data, z signed by the OR and sized by the two-sided P).
Gene sets (user 2026-11-05: full mainland only, no verdict; the six drivers of the gene-based analysis, GENE_DRIVERS,
checked against its leave-one-variant-out table): how many WGS carriers of each driver variant the imputed
arrays recover (dosage >= 0.5), in cases, read against the recovery of case carriers of the rare untyped variants
within 500 kb of the drivers (array_crosscheck's driver windows) with
the same number of control carriers (06.imputed/recovery_by_group.tsv). The controls are in the imputation panel, so
their recovery is listed and carries no information; a typed array call is not used (only the imputed arrays).

Outputs (<out-prefix>.*; every table has a `cohort` column):
  evidence.tsv   one row per cohort x signal x marker x evidence
  signals.tsv    one row per cohort x signal: strongest evidence, verdict (genes: 'carriers counted'), case carriers
                 recovered, every row's verdict
  summary.json   per cohort, the genome-wide numbers and the verdict counts; the imputation; gate G10
  gate.json      gate G10 (array_cohorts.py): full mainland reproduces check/array_crosscheck
"""
IMPUTATION = ('Minimac4 v4.1.4 (VCF header of the imputed set, 2025-10-28) against the AGP3K whole-genome reference '
              'panel, which holds the controls\' own genomes: only the cases\' imputed genotypes are independent of WGS')
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402
import report_lib as RL  # noqa: E402

ASSOC = ('typed', 'imputed', 'imputed (HIBAG)', 'LD proxy')
R2_MIN = 0.8
GENE_COHORT = 'full_mainland'     # gene carriers are counted here only (user 2026-11-05)
# the driver variants of the gene-based analysis (user 2026-11-05: "driven变体要一致"): the leading variant of each
# tiered M/H set and any other whose removal lowers -log10 P by >= 1 in all six tests (integration/gene_signals,
# gene_signals.lovo.tsv); check/array_crosscheck's own list also held STBD1 p.Ser188Arg, whose removal RAISES -log10 P
GENE_DRIVERS = {'chr4:76306990:G:C': 'STBD1 c.220+1G>C', 'chr4:76309559:G:T': 'STBD1 p.Trp212Cys',
                'chr11:59576952:G:A': 'OSBP p.Pro712Ser', 'chr11:59594182:C:T': 'OSBP p.Arg462Gln',
                'chr12:112890858:C:T': 'RPH3A p.Arg544Cys', 'chr9:5073770:G:T': 'JAK2 p.Val617Phe'}
LOVO = K.PROJECT / 'integration' / 'gene_signals' / 'results' / 'gene_signals.lovo.tsv'


def check_drivers():
    """GENE_DRIVERS against the gene-based leave-one-variant-out table: per gene (M/H), the leading variant (the top one
    in each cohort's test with the smaller P, e.g. SKAT-O for JAK2) and any variant lowering -log10 P by >= 1 in all
    six tests; abort if they differ"""
    d = pd.read_csv(LOVO, sep='\t')
    genes = sorted({v.split()[0] for v in GENE_DRIVERS.values()})
    m = d[(d.stratum == 'moderate_high') & d.gene_symbol.isin(genes) & d.dropped.notna()]
    best = m.loc[m.groupby(['gene_symbol', 'cohort']).p_all.transform('min') == m.p_all]
    lead = set(best[best.top.astype(bool)].dropped)
    t = m.groupby('dropped').delta_log10.agg(['min', 'size'])
    found = lead | set(t[(t['min'] >= 1) & (t['size'] == 6)].index)
    if found != set(GENE_DRIVERS):
        sys.exit(f'ABORT: gene drivers differ from {LOVO}: {sorted(found ^ set(GENE_DRIVERS))}')


def z_of(or_, p):
    or_, p = float(or_), float(p)
    if not (np.isfinite(or_) and np.isfinite(p)) or or_ <= 0:
        return np.nan
    return float(np.sign(np.log(or_)) * stats.norm.isf(p / 2))


def assoc_verdict(or_w, p_w, or_a, p_a):
    """the rule of the array check (user 2026-11-05): contradicted = opposite direction or |z WGS - z array| > 3,
    otherwise consistent; the array P is not part of the rule. Returns the verdict, z WGS - z array and both z."""
    zw, za = z_of(or_w, p_w), z_of(or_a, p_a)
    if not np.isfinite(za):
        return 'not checkable', np.nan, zw, za
    dz = zw - za
    if np.sign(za) != np.sign(zw) or abs(dz) > 3:
        return 'contradicted', dz, zw, za
    return 'consistent', dz, zw, za


def conc_typed(arr, s):
    """concordance of a typed array marker with WGS, cases and controls: 1 - discordant / called on both"""
    n = json.loads((arr / '04.compare' / 'summary.json').read_text())['samples']
    out = {}
    for g in ('cases', 'controls'):
        called = n[g] * (1 - s[f'miss_wgs_{g}']) * (1 - s[f'miss_arr_{g}'])
        out[f'conc_{g}'] = float(1 - s[f'disc_{g}'] / called) if called > 0 else np.nan
    return out


def conc_imputed(arr):
    """best-guess concordance with WGS per variant and group (array_cohorts.py step 06): overall, and among the people
    non-reference on WGS or on the array (the carriers' agreement, which the overall share hides for a rare allele)"""
    c = pd.read_csv(arr / '06.imputed' / 'concordance.tsv', sep='\t')
    return {(r.variant, r.group, k): getattr(r, k) for r in c.itertuples() for k in ('concordance', 'nonref_concordance')}


def alt_or(vid, a1, *ors):
    """odds ratios per the ALT allele of the lead (as every other table of the component): array_crosscheck counts its
    A1, often the minor allele (user 2026-11-05: LINC01146, SLC44A2 and ALPI read per REF before)"""
    flip = str(a1) != vid.split(':')[3]
    return [1 / float(o) if flip and pd.notna(o) and float(o) > 0 else o for o in ors]


def gwas_rows(arr, leads, info_r2):
    """one marker per lead: typed, else imputed (R2 >= 0.8), else LD proxy (r2 >= 0.8), else none"""
    cimp = conc_imputed(arr)
    sig = pd.read_csv(arr / '04.compare' / 'signals.tsv', sep='\t', dtype={'lead': str})
    sig = sig[sig.kind == 'common'].set_index('lead')
    imp = pd.read_csv(arr / '06.imputed' / 'leads.tsv', sep='\t', dtype={'lead': str}).set_index('lead')
    imp = imp[~imp.index.duplicated()]          # one lead can sit under two labels (HLA-DQA1 and 'MHC lead')
    sig = sig[~sig.index.duplicated()]
    rows = []
    for unit, vid in leads:
        base = {'signal': unit, 'kind': 'gwas', 'marker': vid}
        s = sig.loc[vid] if vid in sig.index else None
        r2p = float(s.proxy_r2) if s is not None and pd.notna(s.proxy_r2) else np.nan
        ir2 = info_r2.get(vid, np.nan)
        if s is not None and s.proxy == vid:                                      # 1. the lead typed
            ow, oa = alt_or(vid, s.lead_a1, s.sub_wgs_or_proxy, s.sub_arr_or_proxy)
            v, dz, zw_, za_ = assoc_verdict(ow, s.sub_wgs_p_proxy, oa, s.sub_arr_p_proxy)
            rows.append({**base, 'evidence': 'typed', 'array_marker': vid, 'or_allele': 'ALT', 'or_wgs': ow,
                         'p_wgs': s.sub_wgs_p_proxy, 'or_array': oa, 'p_array': s.sub_arr_p_proxy,
                         'dz': dz, 'z_wgs': zw_, 'z_array': za_, **conc_typed(arr, s), 'verdict': v})
        elif vid in imp.index and str(imp.loc[vid, 'in_imputed_set']) == 'True' and ir2 >= R2_MIN:   # 2. imputed
            i = imp.loc[vid]
            ow, oa = alt_or(vid, i.a1, i.or_wgs, i.or_imputed)
            v, dz, zw_, za_ = assoc_verdict(ow, i.p_wgs, oa, i.p_imputed)
            rows.append({**base, 'evidence': 'imputed', 'array_marker': vid, 'info_r2': ir2, 'or_allele': 'ALT',
                         'or_wgs': ow, 'p_wgs': i.p_wgs, 'or_array': oa, 'p_array': i.p_imputed, 'dz': dz, 'z_wgs': zw_, 'z_array': za_,
                         'conc_cases': cimp.get((vid, 'cases', 'concordance'), np.nan),
                         'conc_controls': cimp.get((vid, 'controls', 'concordance'), np.nan),
                         'nonref_conc_cases': cimp.get((vid, 'cases', 'nonref_concordance'), np.nan),
                         'nonref_conc_controls': cimp.get((vid, 'controls', 'nonref_concordance'), np.nan),
                         'verdict': v})
        elif s is not None and np.isfinite(r2p) and r2p >= R2_MIN:               # 3. an LD proxy
            v, dz, zw_, za_ = assoc_verdict(s.sub_wgs_or_proxy, s.sub_wgs_p_proxy, s.sub_arr_or_proxy, s.sub_arr_p_proxy)
            rows.append({**base, 'evidence': 'LD proxy', 'array_marker': s.proxy, 'proxy_r2': r2p,
                         'or_allele': 'proxy A1',
                         'or_wgs': s.sub_wgs_or_proxy, 'p_wgs': s.sub_wgs_p_proxy, 'or_array': s.sub_arr_or_proxy,
                         'p_array': s.sub_arr_p_proxy, 'dz': dz, 'z_wgs': zw_, 'z_array': za_, **conc_typed(arr, s), 'verdict': v})
        else:                                                                     # 4. none
            rows.append({**base, 'evidence': 'none', 'verdict': 'not checkable',
                         'note': 'not typed, not imputed with R2 >= 0.8, no proxy with r2 >= 0.8'})
    return rows


def hla_rows(arr):
    # DQB1:57:V from HIBAG (array_hibag.py, user 2026-11-05); the ridge prediction of array_cohorts.py stays in
    # 04.compare/summary.json for gate G10 only and is not read here
    h = json.loads((arr / '07.hibag' / 'summary.json').read_text())
    t, p = h['assoc_typed'], h['assoc_hibag']
    v, dz, zw_, za_ = assoc_verdict(t['or'], t['p'], p['or'], p['p'])
    r = h['run']
    return [{'signal': RL.V57, 'kind': 'hla', 'marker': RL.V57, 'evidence': 'imputed (HIBAG)', 'or_allele': 'V',
             'array_marker': f"HIBAG Asian model, {r['model_snps_on_array']} of {r['model_snps']} model SNPs on the arrays",
             'or_wgs': t['or'], 'p_wgs': t['p'], 'or_array': p['or'], 'p_array': p['p'], 'dz': dz, 'z_wgs': zw_, 'z_array': za_,
             'case_freq_wgs': h['cases']['v57_freq_typed'], 'case_freq_array': h['cases']['v57_freq_hibag'],
             'control_freq_wgs': h['controls']['v57_freq_typed'], 'control_freq_array': h['controls']['v57_freq_hibag'],
             'r2_cases': h['cases']['v57_r2'], 'r2_controls': h['controls']['v57_r2'],
             'conc_cases': h['cases']['v57_best_agreement'], 'conc_controls': h['controls']['v57_best_agreement'],
             'posterior_median': h['controls']['posterior_median'], 'verdict': v}]


def gene_rows(arr, genes, info_r2):
    """full mainland only: the WGS case carriers of each driver recovered by the imputed arrays (dosage >= 0.5), with the
    recovery expected for a rare variant within 500 kb with as many control carriers; controls listed (panel members)"""
    dr = pd.read_csv(arr / '06.imputed' / 'drivers.tsv', sep='\t')
    dr = dr[dr.id.isin(GENE_DRIVERS)]                    # the gene-based drivers only
    missing = set(GENE_DRIVERS) - set(dr.id)
    if missing:
        sys.exit(f'ABORT: {arr}: no drivers.tsv row for {sorted(missing)}')
    cimp = conc_imputed(arr)
    rows = []
    for var, d in dr.groupby('variant', sort=False):
        base = {'signal': var.split()[0], 'kind': 'gene', 'marker': var}
        ca, co = d[d.group == 'cases'].iloc[0], d[d.group == 'controls'].iloc[0]
        if str(ca.in_imputed_set) != 'True':
            rows.append({**base, 'evidence': 'none', 'verdict': 'not checkable', 'carriers_wgs': int(ca.wgs_carriers),
                         'controls_wgs': int(co.wgs_carriers), 'note': 'not in the imputed arrays'})
            continue
        f, w = int(ca.both), int(ca.wgs_carriers)
        rows.append({**base, 'evidence': 'imputed', 'info_r2': info_r2.get(var, np.nan), 'carriers_wgs': w,
                     'carriers_found': f, 'expected_recovery': float(ca.expected_case_recovery),
                     'carriers_expected': w * float(ca.expected_case_recovery),
                     'mean_dosage': float(ca.mean_dosage_wgs_carriers), 'carriers_extra': int(ca.imputed_only),
                     'controls_wgs': int(co.wgs_carriers), 'controls_found': int(co.both),
                     'conc_cases': cimp.get((ca.id, 'cases', 'concordance'), np.nan),
                     'conc_controls': cimp.get((ca.id, 'controls', 'concordance'), np.nan),
                     'nonref_conc_cases': cimp.get((ca.id, 'cases', 'nonref_concordance'), np.nan),
                     'carriers': f'cases {f} of {w}', 'verdict': 'counted'})
    return [r for r in rows if r['signal'] in genes]


def summarize(ev):
    out = []
    for sig, g in ev.groupby('signal', sort=False):
        best = g.iloc[0]
        vs = g.verdict.tolist()
        rec = ''
        if 'contradicted' in vs:
            overall = 'contradicted'
        elif best.kind != 'gene':
            overall = {'consistent': 'consistent on arrays'}.get(best.verdict, best.verdict)
        else:
            c = g[g.verdict == 'counted']
            overall = 'carriers counted' if len(c) else 'not checkable'
            rec = f'cases {int(c.carriers_found.sum())} of {int(c.carriers_wgs.sum())}' if len(c) else ''
        out.append({'signal': sig, 'kind': best.kind, 'strongest_evidence': best.evidence, 'array_verdict': overall,
                    'case_carriers_recovered': rec,
                    'verdicts': '; '.join(f'{m} {e}: {v}' if best.kind == 'gene' else f'{e}: {v}'
                                          for m, e, v in zip(g.marker, g.evidence, g.verdict))})
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohorts', required=True, help='the per-cohort tables of array_cohorts.py (<dir>/<cohort>/)')
    ap.add_argument('--signals', required=True, help='02.targets/signal_variants.tsv')
    ap.add_argument('--out-prefix', required=True)
    a = ap.parse_args()
    root = Path(a.cohorts)
    check_drivers()
    gate = json.loads((root / 'array.gate.json').read_text())
    if not gate['pass']:
        sys.exit('ABORT: gate G10 did not pass; no array evidence is read')
    Path(f'{a.out_prefix}.gate.json').write_text(json.dumps(gate, indent=2, default=float) + '\n')

    sv = pd.read_csv(a.signals, sep='\t', dtype=str)
    leads = RL.gwas_leads(sv)[['unit', 'variant_id']]
    genes = sorted(sv.loc[sv.kind == 'gene', 'unit'].unique())       # the Tier 1-2 genes
    # the imputation R2 of each marker: one value (the imputation ran once; identical in the three cohorts)
    fm = root / 'full_mainland' / '06.imputed'
    info_r2 = {**dict(zip(pd.read_csv(fm / 'leads.tsv', sep='\t').lead, pd.read_csv(fm / 'leads.tsv', sep='\t').info_r2)),
               **dict(zip(pd.read_csv(fm / 'drivers.tsv', sep='\t').variant, pd.read_csv(fm / 'drivers.tsv', sep='\t').r2))}
    evs, sms, rep = [], [], {'imputation': IMPUTATION, 'gate_G10': gate, 'cohorts': {}, 'r2_min': R2_MIN}
    for c in K.COHORTS:
        arr = root / c
        gn = genes if c == GENE_COHORT else []                      # gene carriers: full mainland only
        ev = pd.DataFrame(gwas_rows(arr, list(leads.itertuples(index=False, name=None)), info_r2) + hla_rows(arr) +
                          gene_rows(arr, gn, info_r2))
        missing = (set(leads.unit) | {RL.V57} | set(gn)) - set(ev.signal)
        if missing:
            sys.exit(f'ABORT: {c}: no array evidence row for {sorted(missing)}')
        sm = summarize(ev)
        evs.append(ev.assign(cohort=c))
        sms.append(sm.assign(cohort=c))
        cs = json.loads((arr / '04.compare' / 'summary.json').read_text())
        gw = cs['genome_wide']
        rep['cohorts'][c] = {
            'people': cs['samples'], 'shared_sites': cs['sites']['shared'], 'qc_sites': cs['sites']['qc_pass'],
            'concordance_cases_T7_vs_ASA24v1': cs['concordance_all_shared_cases']['rate'],
            'concordance_controls_HiSeqX_vs_JSA24v1': cs['concordance_all_shared_controls']['rate'],
            'lambda_wgs': gw['lambda_wgs'], 'lambda_array': gw['lambda_array'], 'r_z_wgs_array': gw['r_z'],
            'slope_z_array_on_wgs': gw['slope_z_array_on_wgs'], 'abs_dz_gt4': gw['abs_dz_gt4'],
            'sites_wgs_p_lt_1e-4': gw['wgs_p_lt_1e-4'],
            'of_which_array_same_direction_p_lt_0.01': gw['wgs_p_lt_1e-4_array_p_lt_1e-2'],
            'signals_by_verdict': sm.array_verdict.value_counts().to_dict()}
        if c == GENE_COHORT:
            g = ev[ev.verdict == 'counted']
            n = ev[(ev.kind == 'gene') & (ev.verdict == 'not checkable')]
            rep['gene_carriers'] = {
                'cohort': c, 'drivers_counted': int(len(g)), 'case_carriers_wgs': int(g.carriers_wgs.sum()),
                'case_carriers_recovered': int(g.carriers_found.sum()),
                'case_carriers_expected': round(float(g.carriers_expected.sum()), 2),
                'case_noncarriers_imputed_as_carriers': int(g.carriers_extra.sum()),
                'not_in_imputed_arrays': n.marker.tolist(),
                'expected': 'recovery of case carriers of the rare untyped variants within 500 kb of the drivers '
                            'with as many control carriers',
                'controls': 'in the imputation panel: their recovery carries no information'}
    rep['coverage'] = ('each cohort\'s unrelated DNBSeq-T7 cases and HiSeq X controls with arrays; NovaSeq and '
                       'DNBSeq-G400RS cases have none')
    first = lambda d: d[['cohort'] + [k for k in d.columns if k != 'cohort']]  # noqa: E731
    ev, sm = first(pd.concat(evs, ignore_index=True)), first(pd.concat(sms, ignore_index=True))
    ev.to_csv(f'{a.out_prefix}.evidence.tsv', sep='\t', index=False)
    sm.to_csv(f'{a.out_prefix}.signals.tsv', sep='\t', index=False)
    Path(f'{a.out_prefix}.summary.json').write_text(json.dumps(rep, indent=2, default=float) + '\n')
    print(sm.pivot(index='signal', columns='cohort', values='array_verdict')[K.COHORTS].to_string())


if __name__ == '__main__':
    main()
