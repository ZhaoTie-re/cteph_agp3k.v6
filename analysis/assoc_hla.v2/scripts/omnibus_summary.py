#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The omnibus tests of residue positions (a main analysis; fixed model, Bonferroni) of every
#           case set and class, stacked and set beside the single-marker scans.
#
#           Each unit (residue position or HLA gene) is compared with the best
#           single marker it contains in the FIXED model (same model, same
#           samples: only the unit of testing differs) and in the RANDOM model
#           (the primary scans: the model differs too), and with the random
#           model's independent single-marker signals (STEPWISE): which unit
#           holds each signal lead.
# Inputs  : omnibus.<class>.tsv, omnibus_signals.<class>.tsv,
#           omnibus_steps.<class>.tsv of every case set (OMNIBUS);
#           <case_set>/_comparison/assoc_all.tsv (COMPARE);
#           _conditional/signals_all.tsv (CONDITIONAL_COMPARE)
# Outputs : omnibus_all.tsv       every unit x cohort, all case sets and classes
#           omnibus_units.tsv     one row per (case set, class, unit): omnibus P
#                                 and calls per cohort, n/3, the best single
#                                 marker per model, the single-marker signal
#                                 leads the unit holds
#           omnibus_signals_all.tsv, omnibus_steps_all.tsv
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import vocab as V                      # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--omnibus', nargs='+', required=True)
    p.add_argument('--signals', nargs='+', required=True)
    p.add_argument('--steps', nargs='+', required=True)
    p.add_argument('--assoc', nargs='+', required=True, help='assoc_all.tsv of every case set')
    p.add_argument('--single-signals', required=True, help='_conditional/signals_all.tsv')
    p.add_argument('--cohorts', required=True)
    return p.parse_args()


def unit_of(marker, cls):
    gene, pos, _ = V.parse_marker(marker, cls)
    return f'{gene}:{pos}' if cls == 'residue' else gene


def main():
    a = parse_args()
    cohorts = a.cohorts.split(',')
    om = pd.concat([pd.read_csv(f, sep='\t', dtype={'unit': str, 'position': str}) for f in a.omnibus],
                   ignore_index=True)
    sg = [pd.read_csv(f, sep='\t', dtype={'lead': str, 'members': str}) for f in a.signals]
    sg = pd.concat([d for d in sg if len(d)], ignore_index=True)
    st = pd.concat([pd.read_csv(f, sep='\t', dtype={'lead': str}) for f in a.steps], ignore_index=True)
    asc = pd.concat([pd.read_csv(f, sep='\t', dtype={'marker': str}) for f in a.assoc], ignore_index=True)
    ss = pd.read_csv(a.single_signals, sep='\t', dtype={'lead': str})
    if set(om.cohort) != set(cohorts):
        sys.exit(f'ERROR: omnibus cohorts {sorted(set(om.cohort))} != {cohorts}')

    asc = asc[asc.test_id == asc.marker].copy()          # one row per test, as in the scans
    asc['unit'] = [unit_of(m, c) for m, c in zip(asc.marker, asc['class'])]
    rows = []
    for (cs, cls, u), g in om.groupby(['case_set', 'class', 'unit'], sort=False):
        g = g.set_index('cohort')
        r = {'case_set': cs, 'class': cls, 'unit': u, 'gene': g.gene.iloc[0]}
        for c in cohorts:
            r[f'df_{c}'] = g.at[c, 'df'] if c in g.index else np.nan
            r[f'p_{c}'] = g.at[c, 'p'] if c in g.index else np.nan
            r[f'sig_{c}'] = bool(g.at[c, 'sig_bonf']) if c in g.index else False
        n = sum(r[f'sig_{c}'] for c in cohorts)
        r['bonf_n'], r['bonf'] = n, f'{n}/{len(cohorts)}'
        for model in ('fixed', 'random'):
            s = asc[(asc.case_set == cs) & (asc['class'] == cls) & (asc.model == model) & (asc.unit == u)
                    & asc.p.notna()]
            if s.empty:
                r[f'best_{model}'] = ''
                continue
            best = s.groupby('marker').p.max().idxmin()      # the marker whose largest P over the cohorts is smallest
            r[f'best_{model}'] = best
            for c in cohorts:
                pc = s[(s.marker == best) & (s.cohort == c)].p
                r[f'best_{model}_p_{c}'] = float(pc.iloc[0]) if len(pc) else np.nan
        leads = ss[(ss.case_set == cs) & (ss['class'] == cls) & (ss.model == 'random')].lead
        r['random_signal_leads'] = ';'.join(m for m in leads if unit_of(m, cls) == u)
        rows.append(r)
    units = pd.DataFrame(rows)

    sg = sg.copy()

    def leads_in(cs, cls, units_):
        return ';'.join(m for m in ss[(ss.case_set == cs) & (ss['class'] == cls) & (ss.model == 'random')].lead
                        if unit_of(m, cls) in units_)
    sg['random_signal_leads_in_unit'] = [leads_in(cs, cls, {u}) for cs, cls, u in zip(sg.case_set, sg['class'], sg.lead)]
    sg['random_signal_leads_in_members'] = [leads_in(cs, cls, set(str(m).split(';')))
                                            for cs, cls, m in zip(sg.case_set, sg['class'], sg.members)]
    om.to_csv('omnibus_all.tsv', sep='\t', index=False, float_format='%.6g')
    units.to_csv('omnibus_units.tsv', sep='\t', index=False, float_format='%.6g')
    sg.to_csv('omnibus_signals_all.tsv', sep='\t', index=False, float_format='%.6g')
    st.to_csv('omnibus_steps_all.tsv', sep='\t', index=False, float_format='%.6g')
    for (cs, cls), g in sg.groupby(['case_set', 'class']):
        print(f'[omnibus_summary] {cs} {cls}: ' + ', '.join(
            f"{r.lead} {r.consistency} (random-model single-marker lead in it: {r.random_signal_leads_in_unit or 'none'}; "
            f"in its members: {r.random_signal_leads_in_members or 'none'})"
            for r in g.itertuples()))


if __name__ == '__main__':
    main()
