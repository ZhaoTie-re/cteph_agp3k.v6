#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The conditional analyses side by side, the random model (the
#           primary model, as in the genome-wide scan) first, every result read
#           across the three cohorts.
#
#   signals_all.tsv          the independent signals of every (case set, model,
#                            class): stepwise.py, selected jointly over the
#                            three cohorts, with their cross-cohort consistency
#                            at entry and marginally, members and carrier alleles
#   steps_all.tsv            every stepwise step
#   hla_given_lead_all.tsv   every HLA marker conditioned on each GWAS lead
#   lead_given_hla_all.tsv   each GWAS lead conditioned on HLA markers
#   ld_lead_hla.tsv          r and r^2 of every HLA marker with each lead
#                            (case set 'all', all samples and PHENO1 = 0)
#   summary.tsv              GWAS lead against the HLA markers, one row per
#                            (case set, cohort, model)
#   consistency.tsv          one row per (case set, model): the signals and
#                            the GWAS-lead results, narrow / intermediate / full
#
#           Tables only: the figures of these analyses are drawn by figures.py
#           (03.conditioning, 04.gwas_link).
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import vocab as V                        # noqa: E402

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--steps', nargs='+', required=True)
    p.add_argument('--signals', nargs='+', required=True)
    p.add_argument('--hla-given', nargs='+', required=True)
    p.add_argument('--lead-given', nargs='+', required=True)
    p.add_argument('--markers', nargs='+', required=True, help="pooled_markers.tsv of every (case set, cohort)")
    p.add_argument('--gwas-leads', required=True)
    p.add_argument('--out-dir', default='.')
    return p.parse_args()


def cat(files, **kw):
    return pd.concat([pd.read_csv(f, sep='\t', **kw) for f in files], ignore_index=True)


def order_frame(d, cols=('case_set', 'cohort', 'model')):
    rank = {**{c: i for i, c in enumerate(V.CASE_SETS)}, **{c: i for i, c in enumerate(V.COHORTS)},
            **{m: i for i, m in enumerate(('random', 'fixed'))}, **{c: i for i, c in enumerate(V.CLASSES)}}
    return d.sort_values(list(cols), key=lambda s: s.map(rank)).reset_index(drop=True)


def main():
    a = parse_args()
    out = Path(a.out_dir)
    steps = order_frame(cat(a.steps, dtype={'lead': str}), cols=('case_set', 'model', 'class'))
    sig = cat(a.signals, dtype={'lead': str})
    sig = order_frame(sig, cols=('case_set', 'model', 'class')) if len(sig) else sig
    hg = cat(a.hla_given, dtype={'marker': str})
    lg = cat(a.lead_given, dtype={'conditioned_on': str, 'position': str, 'position_reference': str})
    mk = cat(a.markers, dtype={'marker': str})
    gl = pd.read_csv(a.gwas_leads, sep='\t')
    leads = gl.loc[gl.status == 'kept', 'lead_id'].tolist()

    steps.to_csv(out / 'steps_all.tsv', sep='\t', index=False, float_format='%.6g')
    sig.to_csv(out / 'signals_all.tsv', sep='\t', index=False, float_format='%.6g')
    first = {}                                     # (case set, model, class) -> first signal's lead
    for (cs, mod, c), d in sig.groupby(['case_set', 'model', 'class']):
        first[(cs, mod, c)] = d.sort_values('signal').lead.iloc[0]
    order_frame(hg).to_csv(out / 'hla_given_lead_all.tsv', sep='\t', index=False, float_format='%.6g')
    order_frame(lg).to_csv(out / 'lead_given_hla_all.tsv', sep='\t', index=False, float_format='%.6g')
    ld = mk[(mk.case_set == 'all') & (mk['class'] != 'gwas_lead')].copy()
    keep = ['cohort', 'marker', 'class', 'gene', 'position', 'residue']
    ldl = []
    for ldid in leads:
        t = ld[keep + [f'r_{ldid}', f'r_group0_{ldid}']].rename(
            columns={f'r_{ldid}': 'r', f'r_group0_{ldid}': 'r_controls'})
        t.insert(0, 'lead', ldid)
        t['r2'] = t['r'] ** 2
        ldl.append(t)
    ldt = pd.concat(ldl, ignore_index=True)
    ldt['position'] = ldt['position'].astype('Int64')
    ldt.sort_values(['lead', 'cohort', 'r2'], ascending=[True, True, False]).to_csv(
        out / 'ld_lead_hla.tsv', sep='\t', index=False, float_format='%.6g')

    # ── GWAS lead against the HLA markers, per (case set, cohort, model) ────
    rows = []
    for (cs, coh, mod), d0 in lg.groupby(['case_set', 'cohort', 'model'], sort=False):
        rec = {'case_set': cs, 'cohort': coh, 'model': mod,
               'residue_signal1': first.get((cs, mod, 'residue'), ''),
               'allele_signal1': first.get((cs, mod, 'allele'), '')}
        for ldid in leads:
            d = d0[d0.lead == ldid]
            single = d[~d.why.str.contains('jointly')]
            joint = d[d.why == 'signal_leads_jointly']
            pos = d[d.why.str.startswith('position_jointly')]
            p1 = rec['residue_signal1'].rsplit(':', 1)[0] if rec['residue_signal1'] else ''
            g1 = pos[pos.position == p1]
            rec[f'{ldid}|position1'] = p1
            rec[f'{ldid}|p_given_position1'] = g1.lead_p_given.iloc[0] if len(g1) else np.nan
            go = pos[pos.why.str.contains('omnibus_lead')]
            rec[f'{ldid}|omnibus_position1'] = go.position.iloc[0] if len(go) else ''
            rec[f'{ldid}|p_given_omnibus_position1'] = go.lead_p_given.iloc[0] if len(go) else np.nan
            rec[f'{ldid}|marginal_p'] = d.lead_p_marginal.iloc[0]
            for c in V.CLASSES:
                g = single[single.conditioned_on == rec[f'{c}_signal1']]
                rec[f'{ldid}|p_given_{c}_signal1'] = g.lead_p_given.iloc[0] if len(g) else np.nan
            rec[f'{ldid}|p_given_all_signals'] = joint.lead_p_given.iloc[0] if len(joint) else np.nan
            b = single.loc[single.lead_p_given.idxmax()]
            rec[f'{ldid}|most_explaining'] = b.conditioned_on
            rec[f'{ldid}|p_given_most_explaining'] = b.lead_p_given
            t = single.loc[single.r.abs().idxmax()]
            rec[f'{ldid}|best_ld_tag'] = t.conditioned_on
            rec[f'{ldid}|best_ld_tag_r2'] = t.r ** 2
            rec[f'{ldid}|p_given_best_ld_tag'] = t.lead_p_given
            h = hg[(hg.case_set == cs) & (hg.cohort == coh) & (hg.model == mod) & (hg.lead == ldid)]
            rec[f'{ldid}|hla_min_p_given_lead'] = h.p_given_lead.min()
            rec[f'{ldid}|hla_top_given_lead'] = h.loc[h.p_given_lead.idxmin(), 'marker'] if h.p_given_lead.notna().any() else ''
            r1 = h[h.marker == rec['residue_signal1']]
            rec[f'{ldid}|residue_signal1_p_given_lead'] = r1.p_given_lead.iloc[0] if len(r1) else np.nan
        rows.append(rec)
    summ = order_frame(pd.DataFrame(rows))
    summ.to_csv(out / 'summary.tsv', sep='\t', index=False, float_format='%.6g')

    # ── per (case set, model), narrow / intermediate / full ─────────────────
    fmt = lambda xs: '/'.join('-' if pd.isna(x) else f'{x:.1e}' for x in xs)        # noqa: E731
    crow = []
    for (cs, mod), d in summ.groupby(['case_set', 'model'], sort=False):
        d = d.set_index('cohort').reindex(list(V.COHORTS))
        rec = {'case_set': cs, 'model': mod}
        for c in V.CLASSES:
            e = sig[(sig.case_set == cs) & (sig.model == mod) & (sig['class'] == c)] if len(sig) else sig
            rec[f'{c}_signals'] = '; '.join(f'{r.lead} ({r.gec})' for r in e.itertuples()) or '-'
        for ldid in leads:
            rec[f'{ldid}|marginal_p<1e-5_all'] = bool((d[f'{ldid}|marginal_p'] < 1e-5).all())
            rec[f'{ldid}|marginal_p'] = fmt(d[f'{ldid}|marginal_p'])
            rec[f'{ldid}|p_given_residue_signal1'] = fmt(d[f'{ldid}|p_given_residue_signal1'])
            rec[f'{ldid}|residue_signal1_p_given_lead'] = fmt(d[f'{ldid}|residue_signal1_p_given_lead'])
            rec[f'{ldid}|most_explaining'] = '/'.join(d[f'{ldid}|most_explaining'].fillna('-'))
            rec[f'{ldid}|p_given_most_explaining'] = fmt(d[f'{ldid}|p_given_most_explaining'])
            rec[f'{ldid}|best_ld_tag'] = '/'.join(d[f'{ldid}|best_ld_tag'].fillna('-'))
            rec[f'{ldid}|p_given_best_ld_tag'] = fmt(d[f'{ldid}|p_given_best_ld_tag'])
        crow.append(rec)
    cons = pd.DataFrame(crow)
    cons.to_csv(out / 'consistency.tsv', sep='\t', index=False)

    print(f'[conditional_compare] {len(sig)} signals over {steps.groupby(["case_set", "model", "class"]).ngroups} '
          f'(case set, model, class) walks; {len(summ)} GWAS-lead analyses')


if __name__ == '__main__':
    main()
