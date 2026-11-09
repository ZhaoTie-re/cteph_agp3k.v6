#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Every HLA marker (alleles and residues) conditioned on each
#           independent-signal lead of the case set, and on each whole position
#           that conditional.py conditions the GWAS lead on, for one (case set,
#           cohort, model), with the SAME engine as the marginal scan (engine.py)
#           on the pooled PLINK set. It completes the HLA | signal side that the
#           stepwise selection leaves open: STEPWISE only steps within a class
#           that has a signal, so e.g. the alleles are never retested given the
#           residue lead. The GWAS lead itself is reported too.
#
#           Conditioning sets: each signal lead alone (both classes) and each
#           whole position (every residue column of the position but its most
#           frequent; conditional.position_sets, so the sets are identical to
#           the ones the GWAS lead is conditioned on).
# Inputs  : as conditional.py
# Outputs : hla_given_signal.tsv   one row per (conditioning set, marker)
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from conditional import parse_args, position_sets      # noqa: E402
from engine import Engine                                # noqa: E402


def main():
    a = parse_args()
    out = Path(a.out_dir)
    mk = pd.read_csv(a.markers, sep='\t', dtype={'marker': str})
    coding = json.load(open(a.coding))
    order = dict(zip(mk['marker'], mk['file_order']))
    hla = mk.loc[mk['class'] != 'gwas_lead', 'marker'].tolist()
    leads = mk.loc[mk['class'] == 'gwas_lead', 'marker'].tolist()
    cls = dict(zip(mk['marker'], mk['class']))
    eng = Engine(a, order, coding[a.model], prefix='sig')
    key = {'case_set': a.case_set, 'cohort': a.cohort, 'model': a.model}
    signal_leads = [ln.strip() for f in a.signal_leads for ln in open(f) if ln.strip()]
    signal_leads = [m for m in signal_leads if m in order]        # an allele may be untested in this cohort
    positions = position_sets(a, hla, cls, signal_leads)
    res0 = eng.run([])

    sets = [([m], 'signal_lead', '') for m in signal_leads]
    sets += [(cset, w, u) for u, (cset, w, _ref) in positions.items()]
    rows = []
    for cset, why, unit in sets:
        rc = eng.run(cset)
        one = cset[0] if len(cset) == 1 else ';'.join(cset)
        for m_ in hla + leads:
            rows.append(dict(key, conditioned_on=one, why=why, position=unit, marker=m_, **{'class': cls[m_]},
                             p_marginal=res0.at[m_, 'p'] if m_ in res0.index else np.nan,
                             or_marginal=res0.at[m_, 'or'] if m_ in res0.index else np.nan,
                             p_given=rc.at[m_, 'p'] if m_ in rc.index else np.nan,
                             or_given=rc.at[m_, 'or'] if m_ in rc.index else np.nan,
                             errcode=rc.at[m_, 'errcode'] if m_ in rc.index else 'NO_RESULT'))
    cols = ['case_set', 'cohort', 'model', 'conditioned_on', 'why', 'position', 'marker', 'class', 'p_marginal',
            'or_marginal', 'p_given', 'or_given', 'errcode']
    pd.DataFrame(rows, columns=cols).to_csv(out / 'hla_given_signal.tsv', sep='\t', index=False, float_format='%.6g')
    print(f"[hla_given_signal] {a.case_set} {a.cohort} {a.model}: {len(sets)} conditioning set(s) "
          f"({', '.join(s[1] + ':' + (s[2] or s[0][0]) for s in sets) or 'none'}); {eng.n_runs} engine runs")


if __name__ == '__main__':
    main()
