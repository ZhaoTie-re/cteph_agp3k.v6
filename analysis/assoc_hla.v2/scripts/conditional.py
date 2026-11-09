#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The relation of the genome-wide scan's MHC lead(s) to the HLA
#           markers, for one (case set, cohort, model), with the SAME engine as
#           the marginal scan (engine.py) on the pooled PLINK set (alleles +
#           residues + the GWAS MHC lead[s]; build_pooled.py).
#
#   1  HLA | LEAD: every HLA marker conditioned on each GWAS lead.
#   2  LEAD | HLA: each GWAS lead conditioned on one HLA marker at a time, for
#      the markers most correlated with it (top --n-tag by |r|), the markers
#      with the smallest marginal P (top --n-tag) and the independent-signal
#      leads of both classes (stepwise.py, same case set and model); on all
#      signal leads jointly; and on WHOLE POSITIONS: every residue column of a
#      position but its most frequent (a position's copies sum to 2), for the
#      positions of the residue signal leads and the omnibus position leads
#      (OMNIBUS, same case set). The reference residue is the most frequent in
#      the pooled set, so both models condition on the same columns.
#
#           The independent HLA signals themselves are selected jointly over
#           the three cohorts by stepwise.py; this script only reads their
#           leads. P, collinearity and the gates are engine.py's.
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from engine import Engine              # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case-set', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--model', required=True, choices=('fixed', 'random'))
    p.add_argument('--bfile', required=True, help='pooled PLINK prefix')
    p.add_argument('--markers', required=True, help='pooled_markers.tsv')
    p.add_argument('--coding', required=True, help='pheno_coding.json')
    p.add_argument('--pheno-fixed', help='fixed model: pheno_covar.fixed.tsv')
    p.add_argument('--keep', help='fixed model: fixed.keep')
    p.add_argument('--covar-name', default='SEX,PC1_AVG-PC10_AVG')
    p.add_argument('--rda', help='random model: null.rda')
    p.add_argument('--vr', help='random model: null.varianceRatio.txt')
    p.add_argument('--plink2', required=True)
    p.add_argument('--saige-rscript', required=True)
    p.add_argument('--saige-step2', required=True)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--signal-leads', nargs='*', default=[], help='leads.<class>.txt from stepwise.py')
    p.add_argument('--position-leads', nargs='*', default=[],
                   help='omnibus_signals.residue.tsv from OMNIBUS: the omnibus position leads')
    p.add_argument('--n-tag', type=int, default=10)
    p.add_argument('--out-dir', default='.')
    return p.parse_args()


def residue_freqs(a, markers):
    """Copy frequency of each residue marker (the ALT, 'present') over the pooled set's samples."""
    Path('freq.ids').write_text('\n'.join(markers) + '\n')
    subprocess.run([a.plink2, '--bfile', a.bfile, '--extract', 'freq.ids', '--freq', '--out', 'freq'], check=True,
                   stdout=subprocess.DEVNULL)
    f = pd.read_csv('freq.afreq', sep='\t', dtype={'ID': str})
    return dict(zip(f['ID'], f['ALT_FREQS']))


def position_sets(a, hla, cls, signal_leads):
    """{position: (columns, why, reference)} for the residue signal leads' positions and the omnibus leads."""
    src = {}
    for m_ in signal_leads:
        if cls.get(m_) == 'residue':
            src.setdefault(m_.rsplit(':', 1)[0], []).append('residue_lead')
    for f in a.position_leads:
        for u in pd.read_csv(f, sep='\t', dtype={'lead': str}).lead.dropna():
            src.setdefault(u, []).append('omnibus_lead')
    cols = {u: [m_ for m_ in hla if cls[m_] == 'residue' and m_.rsplit(':', 1)[0] == u] for u in src}
    cols = {u: c for u, c in cols.items() if c}          # a position may be untested in this cohort
    freq = residue_freqs(a, sorted({m_ for c in cols.values() for m_ in c})) if cols else {}
    out = {}
    for u, c in cols.items():
        ref = max(c, key=lambda m_: freq[m_]) if len(c) > 1 else ''
        out[u] = ([m_ for m_ in c if m_ != ref], 'position_jointly:' + ','.join(dict.fromkeys(src[u])), ref)
    return out


def main():
    a = parse_args()
    out = Path(a.out_dir)
    mk = pd.read_csv(a.markers, sep='\t', dtype={'marker': str})
    coding = json.load(open(a.coding))
    order = dict(zip(mk['marker'], mk['file_order']))
    hla = mk.loc[mk['class'] != 'gwas_lead', 'marker'].tolist()
    leads = mk.loc[mk['class'] == 'gwas_lead', 'marker'].tolist()
    cls = dict(zip(mk['marker'], mk['class']))
    eng = Engine(a, order, coding[a.model])
    key = {'case_set': a.case_set, 'cohort': a.cohort, 'model': a.model}
    signal_leads = [ln.strip() for f in a.signal_leads for ln in open(f) if ln.strip()]
    signal_leads = [m for m in signal_leads if m in order]        # an allele may be untested in this cohort
    positions = position_sets(a, hla, cls, signal_leads)
    res0 = eng.run([])

    # ── 1. HLA | lead ───────────────────────────────────────────────────────
    hg_rows = []
    given = {}
    for ld in leads:
        rg = eng.run([ld])
        given[ld] = rg
        for m_ in hla:
            hg_rows.append(dict(key, lead=ld, marker=m_, **{'class': cls[m_]},
                                r=float(mk.set_index('marker').at[m_, f'r_{ld}']),
                                p_marginal=res0.at[m_, 'p'] if m_ in res0.index else np.nan,
                                or_marginal=res0.at[m_, 'or'] if m_ in res0.index else np.nan,
                                p_given_lead=rg.at[m_, 'p'] if m_ in rg.index else np.nan,
                                or_given_lead=rg.at[m_, 'or'] if m_ in rg.index else np.nan,
                                errcode=rg.at[m_, 'errcode'] if m_ in rg.index else 'NO_RESULT'))
    hg = pd.DataFrame(hg_rows)
    hg.to_csv(out / 'hla_given_lead.tsv', sep='\t', index=False, float_format='%.6g')

    # ── 2. lead | HLA ───────────────────────────────────────────────────────
    lg_rows = []
    m_idx = mk.set_index('marker')
    for ld in leads:
        rr = m_idx.loc[hla, f'r_{ld}'].abs().dropna().sort_values(ascending=False)
        pp = res0.loc[res0.index.intersection(hla), 'p'].dropna().sort_values()
        why = {}
        for m_ in rr.index[:a.n_tag]:
            why.setdefault(m_, []).append('ld_tag')
        for m_ in pp.index[:a.n_tag]:
            why.setdefault(m_, []).append('top_p')
        for m_ in signal_leads:
            why.setdefault(m_, []).append('signal_lead')
        sets = [([m_], ','.join(w), '', '') for m_, w in why.items()]
        if signal_leads:
            sets.append((signal_leads, 'signal_leads_jointly', '', ''))
        sets += [(cset, w, u, ref) for u, (cset, w, ref) in positions.items()]
        for cset, w, u, ref in sets:
            rc = eng.run(cset)
            one = cset[0] if len(cset) == 1 else ';'.join(cset)
            lg_rows.append(dict(key, lead=ld, conditioned_on=one, why=w, position=u, position_reference=ref,
                                r=float(m_idx.at[cset[0], f'r_{ld}']) if len(cset) == 1 else np.nan,
                                lead_p_marginal=float(res0.at[ld, 'p']), lead_or_marginal=float(res0.at[ld, 'or']),
                                lead_p_given=float(rc.at[ld, 'p']), lead_or_given=float(rc.at[ld, 'or']),
                                marker_p_marginal=float(res0.at[cset[0], 'p']) if len(cset) == 1 else np.nan,
                                marker_p_given_lead=(float(given[ld].at[cset[0], 'p'])
                                                     if len(cset) == 1 and cset[0] in given[ld].index else np.nan)))
    lg = pd.DataFrame(lg_rows)
    lg.to_csv(out / 'lead_given_hla.tsv', sep='\t', index=False, float_format='%.6g')

    summ = dict(key, n_hla=len(hla), n_engine_runs=eng.n_runs, signal_leads=signal_leads, leads=leads,
                positions={u: {'columns': c, 'why': w, 'reference': r} for u, (c, w, r) in positions.items()})
    json.dump(summ, open(out / 'conditional.json', 'w'), indent=1, default=str)
    print(f"[conditional] {a.case_set} {a.cohort} {a.model}: signal leads {signal_leads}; positions "
          f"{list(positions)}; {eng.n_runs} engine runs")
    for ld in leads:
        d = lg[(lg.lead == ld) & (lg.position == '')]
        best = d.loc[d.lead_p_given.idxmax()]
        print(f"[conditional]   {ld}: marginal P {d.lead_p_marginal.iloc[0]:.2e}; largest P given one "
              f"HLA set {best.lead_p_given:.2e} ({best.conditioned_on}, r {best.r:.2f})")


if __name__ == '__main__':
    main()
