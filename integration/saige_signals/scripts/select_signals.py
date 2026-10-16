#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Decide what the table shows. Loci = the assoc_saige cross-cohort
#           panels (b: genome-wide in >= 1 cohort; c: suggestive in every
#           cohort and genome-wide in none), each with its representative lead,
#           any other cohort lead of the same locus, and the union over cohorts
#           of its SuSiE credible-set members. Every locus x cohort pair is
#           mapped to the peak that cohort called there (or to none).
# Inputs  : lead_crosscohort.tsv, cs_variants_all.tsv, <cohort>/03.peaks/peaks.tsv
# Outputs : loci.tsv        one row per locus (order, label, class)
#           variants.tsv    one row per table variant
#           sel.tsv         the selected lead rows of lead_crosscohort (reference values)
#           cs_members.tsv  credible-set member rows of the selected peaks
#           pairs.tsv       locus x cohort -> peak (or ''), and the locus's variant list
#           selection.json  counts
# All checks are invariants of the data, never fixed counts.
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sig_common import csv_list, kv_paths, log, num, read_tsv, split_id  # noqa: E402


def load_selection(cc_path, peaks, cohorts):
    cc = read_tsv(cc_path)
    cc['P_num'] = num(cc['P'])
    sel = cc[cc.panel.isin(['b', 'c'])].copy()
    assert len(sel), 'no locus in the cross-cohort panels'
    assert set(sel.cohort) <= set(cohorts), set(sel.cohort) - set(cohorts)
    sel = sel[sel.cohort.isin(cohorts)]
    # every selected variant is reported once per cohort
    assert len(sel) == sel.variant_id.nunique() * len(cohorts), (len(sel), sel.variant_id.nunique())
    assert not sel.duplicated(['variant_id', 'cohort']).any()
    rep = sel[sel.locus_rep == '1']
    assert (rep.groupby('locus_id').variant_id.nunique() == 1).all(), 'a locus has two representative leads'
    assert set(rep.locus_id) == set(sel.locus_id), 'a selected locus has no representative lead'
    alt = set(sel[sel.locus_rep == '0'].variant_id)
    assert not (alt & set(rep.variant_id))

    loci = (rep.groupby('locus_id')
            .agg(variant_id=('variant_id', 'first'), label=('Gene', 'first'),
                 panel=('panel', 'first'), best_tier=('best_tier', 'first'),
                 min_p=('P_num', 'min'))
            .reset_index())
    # the cross-cohort figure's order: panel b first, then the smallest P any cohort reports
    loci = loci.sort_values(['panel', 'min_p'], ascending=[True, True]).reset_index(drop=True)
    loci['order'] = np.arange(1, len(loci) + 1)
    loci['tier'] = loci.best_tier.map({'genome_wide': 'GW', 'suggestive': 'SG'})
    assert loci.tier.notna().all(), set(loci.best_tier)

    # (locus, cohort) -> the peak called there, from EVERY lead of the locus
    # (a cohort may have called its peak under a lead that is not in the panel)
    tiers = {}
    for c in cohorts:
        pk = read_tsv(peaks[c])
        tiers[c] = dict(zip(pk.peak_id, pk.tier))
    all_loc = cc[cc.locus_id.isin(loci.locus_id) & cc.cohort.isin(cohorts)]
    peak_map = {}
    for _, r in all_loc[all_loc.peak_id != '.'].iterrows():
        key = (r.locus_id, r.cohort)
        assert key not in peak_map or peak_map[key][0] == r.peak_id, key
        peak_map[key] = (r.peak_id, tiers[r.cohort][r.peak_id])
    assert {l for l, _c in peak_map} == set(loci.locus_id), 'a locus has no peak in any cohort'
    inv = [(c, p) for (_l, c), (p, _t) in peak_map.items()]
    assert len(inv) == len(set(inv)), 'one peak is mapped to two loci'
    return cc, sel, loci, peak_map


def load_cs_members(cs_path, peak_map):
    cs = read_tsv(cs_path)
    keys = {(c, p) for (_l, c), (p, _t) in peak_map.items()}
    cs = cs[[(c, p) in keys for c, p in zip(cs.cohort, cs.peak_id)]].copy()
    inv = {(p, c): l for (l, c), (p, _t) in peak_map.items()}
    cs['locus_id'] = [inv[(p, c)] for p, c in zip(cs.peak_id, cs.cohort)]
    cs['pip_num'] = num(cs.pip)
    cs['cs_index'] = num(cs.cs).astype(int)         # set number within the cohort's peak (1..k)
    assert (cs.cs_index >= 1).all()
    size = cs.groupby(['cohort', 'peak_id', 'cs_index']).variant_id.transform('size')
    assert (size == num(cs.cs_size)).all(), 'credible-set rows disagree with the recorded set size'
    n_lead_sets = int(cs.groupby(['cohort', 'peak_id']).is_lead.apply(lambda s: (s == 'True').sum()).sum())
    return cs, n_lead_sets


def build_variants(sel, loci, cs, cohorts, codes):
    lead_rows = sel.drop_duplicates('variant_id')[['variant_id', 'locus_id', 'locus_rep', 'rsID',
                                                   'rsID_unresolved', 'Gene', 'EA', 'OA']].copy()
    lead_rows['role'] = np.where(lead_rows.locus_rep == '1', 'Lead', 'Alt lead')
    # annotation must be identical across the cohort rows of a lead
    for col in ['rsID', 'rsID_unresolved', 'Gene', 'EA', 'OA']:
        assert (sel.groupby('variant_id')[col].nunique() == 1).all(), col

    members = cs[~cs.variant_id.isin(lead_rows.variant_id)]
    for col in ['rsID', 'rsID_unresolved', 'Gene', 'EA', 'OA', 'locus_id']:
        assert (members.groupby('variant_id')[col].nunique() <= 1).all(), col
    mem_rows = members.drop_duplicates('variant_id')[['variant_id', 'locus_id', 'rsID', 'rsID_unresolved',
                                                     'Gene', 'EA', 'OA']].copy()
    mem_rows['locus_rep'] = '.'
    mem_rows['role'] = 'CS member'

    v = pd.concat([lead_rows, mem_rows], ignore_index=True)
    assert v.variant_id.is_unique
    roles = v.role.value_counts()
    assert roles.get('Lead', 0) == len(loci)
    assert roles.get('Alt lead', 0) == int((sel.drop_duplicates('variant_id').locus_rep == '0').sum())
    assert roles.get('CS member', 0) == len(set(cs.variant_id) - set(lead_rows.variant_id))

    pipmax = cs.groupby('variant_id').pip_num.max()
    code = dict(zip(cohorts, codes))
    cs_coh = cs.groupby('variant_id').cohort.apply(lambda s: ''.join(code[c] for c in cohorts if c in set(s)))
    v['max_pip'] = v.variant_id.map(pipmax).fillna(-1.0)
    v['cs_cohorts'] = v.variant_id.map(cs_coh).fillna('')
    v = v.merge(loci[['locus_id', 'order', 'label', 'tier', 'panel']], on='locus_id', how='left')
    assert v.order.notna().all()

    parts = v.variant_id.map(split_id)
    v['chrom'] = [p[0] for p in parts]
    v['pos'] = [p[1] for p in parts]
    v['ref'] = [p[2] for p in parts]
    v['alt'] = [p[3] for p in parts]
    assert (v.EA == v.alt).all() and (v.OA == v.ref).all(), 'effect allele is not ALT somewhere'

    role_rank = {'Lead': 0, 'Alt lead': 1, 'CS member': 2}
    v['role_rank'] = v.role.map(role_rank)
    v = v.sort_values(['order', 'role_rank', 'max_pip', 'pos'],
                      ascending=[True, True, False, True]).reset_index(drop=True)
    return v


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--lead-crosscohort', required=True)
    ap.add_argument('--cs-variants', required=True)
    ap.add_argument('--peaks', action='append', required=True, help='cohort=<cohort>/03.peaks/peaks.tsv')
    ap.add_argument('--cohorts', required=True, help='comma list, upstream order')
    ap.add_argument('--cohort-codes', required=True, help='comma list, one letter per cohort')
    ap.add_argument('--out-loci', required=True)
    ap.add_argument('--out-variants', required=True)
    ap.add_argument('--out-sel', required=True)
    ap.add_argument('--out-cs', required=True)
    ap.add_argument('--out-pairs', required=True)
    ap.add_argument('--out-json', required=True)
    a = ap.parse_args()
    cohorts, codes = csv_list(a.cohorts), csv_list(a.cohort_codes)
    assert len(codes) == len(cohorts)
    peaks = kv_paths(a.peaks)
    assert set(peaks) == set(cohorts), (set(peaks), cohorts)

    cc, sel, loci, peak_map = load_selection(a.lead_crosscohort, peaks, cohorts)
    cs, n_lead_sets = load_cs_members(a.cs_variants, peak_map)
    v = build_variants(sel, loci, cs, cohorts, codes)

    pairs = []
    for _, l in loci.iterrows():
        ids = ','.join(sorted(v[v.locus_id == l.locus_id].variant_id))
        for c in cohorts:
            peak, tier = peak_map.get((l.locus_id, c), ('', ''))
            pairs.append({'locus_id': l.locus_id, 'label': l.label, 'cohort': c,
                          'peak_id': peak, 'peak_tier': tier, 'variants': ids})
    pairs = pd.DataFrame(pairs)

    loci.to_csv(a.out_loci, sep='\t', index=False)
    v.to_csv(a.out_variants, sep='\t', index=False)
    sel.drop(columns=['P_num']).to_csv(a.out_sel, sep='\t', index=False)
    cs.to_csv(a.out_cs, sep='\t', index=False)
    pairs.to_csv(a.out_pairs, sep='\t', index=False)
    counts = {
        'n_cohorts': len(cohorts), 'n_loci': len(loci),
        'n_loci_by_class': loci.tier.value_counts().to_dict(),
        'loci_in_order': list(loci.label),
        'n_variants': len(v), 'roles': v.role.value_counts().to_dict(),
        'alt_leads': sorted(v[v.role == 'Alt lead'].variant_id),
        'n_sel_rows': len(sel), 'n_pairs': len(pairs), 'n_pairs_with_peak': int((pairs.peak_id != '').sum()),
        'pairs_without_peak': [f'{r.label}:{r.cohort}' for _, r in pairs[pairs.peak_id == ''].iterrows()],
        'n_cs_rows': len(cs), 'n_cs_variants': int(cs.variant_id.nunique()),
        'n_credible_sets': int(cs.groupby(['cohort', 'peak_id', 'cs_index']).ngroups),
        'n_peaks_with_cs': int(cs.groupby(['cohort', 'peak_id']).ngroups),
        'n_peaks_lead_in_a_set': n_lead_sets,
    }
    with open(a.out_json, 'w') as fh:
        json.dump(counts, fh, indent=2)
    log(f'selection: {len(loci)} loci ({", ".join(loci.label)}); {len(v)} variants '
        f'({counts["roles"]}); {counts["n_pairs_with_peak"]}/{len(pairs)} locus x cohort pairs with a peak')


if __name__ == '__main__':
    main()
