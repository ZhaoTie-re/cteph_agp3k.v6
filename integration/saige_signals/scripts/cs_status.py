#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Credible-set status of every variant x cohort, with the reason when
#           the variant is not in that cohort's set. Combines assoc.tsv (which
#           rows are in a set, the set sizes and cuts) with the SuSiE refits of
#           susie_cs_diagnostic.R (whether a cohort without a set had no effect
#           fitted, or had only impure sets).
#             in_cs            member of set k
#             below_cut        the cohort has a set; PIP below its smallest member
#             no_cs_null       no effect fitted (all PIP = 0)
#             no_cs_impure     effect(s) fitted; every 95 % set had min |r| below threshold
#             not_fine_mapped  no peak in the cohort, or outside its window
#             not_tested       absent from the cohort's scan
# Inputs  : the refits' effects / variants / fit.json files (one per locus x
#           cohort with a peak), assoc.tsv, pairs.tsv, cs_members.tsv, the
#           upstream assoc_saige run manifest
# Outputs : assoc_status.tsv, susie_sets.tsv (one row per fit x effect), cs_status.json
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sig_common import as_bool, load_json, log, num, read_table, read_tsv, sig4, write_tsv  # noqa: E402

EFF_NUM = ('n_gwas', 'n_variants', 'L', 'niter', 'n_effects', 'n_sets_pre', 'n_sets_pure',
           'n_credible_sets_stored', 'max_abs_dpip', 'effect', 'V', 'dedup_of', 'set_size', 'set_coverage',
           'min_abs_corr', 'mean_abs_corr', 'median_abs_corr', 'top_alpha', 'n_table_variants')
EFF_BOOL = ('active', 'purity_pass', 'stored', 'contains_lead', 'converged')


def load_fits(effects, fvars, fjson, pairs, manifest, tol):
    eff = pd.concat([read_tsv(p) for p in effects], ignore_index=True)
    var = pd.concat([read_tsv(p) for p in fvars], ignore_index=True)
    meta = [load_json(p) for p in fjson]
    for c in EFF_NUM:
        eff[c] = num(eff[c])
    for c in EFF_BOOL:
        eff[c] = eff[c].map({'TRUE': True, 'FALSE': False})
    for c in ('pip_stored', 'pip_new', 'alpha_pure', 'alpha_cut_pure', 'nearest_set', 'stored_cs'):
        var[c] = num(var[c])

    # one fit per locus x cohort with a peak, none missing, none extra
    want = {(r.cohort, r.peak_id) for _, r in pairs[pairs.peak_id != ''].iterrows()}
    got = {(m['cohort'], m['locus_id']) for m in meta}
    assert got == want and len(meta) == len(want), (want ^ got)
    fits = eff.drop_duplicates(['cohort', 'locus_id'])
    assert len(fits) == len(want) and len(eff) == int(fits.L.sum())
    s = manifest['susie']
    for m in meta:
        assert m['L'] == s['L'] and abs(m['coverage'] - s['coverage']) < 1e-12 \
            and abs(m['min_abs_corr'] - s['min_abs_corr']) < 1e-12, (m['cohort'], m['locus_id'])
    assert fits.converged.all() and fits.max_abs_dpip.max() < tol, fits.max_abs_dpip.max()
    log(f'SuSiE refits: {len(fits)} fits reproduce the stored PIPs (max |dPIP| = {fits.max_abs_dpip.max():.1e}); '
        f'{int((fits.n_sets_pure >= 1).sum())} with a set, {int((fits.n_effects == 0).sum())} null fits, '
        f'{int(((fits.n_effects > 0) & (fits.n_sets_pure == 0)).sum())} impure')
    return eff, var, meta


def derive_cs_status(assoc, eff, var):
    """Per variant x cohort: in a set, below the cut, or why there is no set."""
    fits = eff.drop_duplicates(['cohort', 'locus_id']).set_index(['cohort', 'locus_id'])
    # the strongest active effect of each fit (largest prior variance V): the dropped 95% set shown for an
    # impure fit (a single-effect fit's only effect; in the multi-effect fits it dominates the near-null rest)
    top = (eff[eff.active].sort_values('V', ascending=False, kind='stable')
           .drop_duplicates(['cohort', 'locus_id']).set_index(['cohort', 'locus_id']))
    vkey = var.set_index(['cohort', 'locus_id', 'ID'])
    out = []
    for _, r in assoc.iterrows():
        rec = {'n_effects': np.nan, 'n_sets_pre': np.nan, 'n_sets_pure': np.nan,
               'dropped_set_size': np.nan, 'dropped_set_min_r': np.nan,
               'alpha_pure': np.nan, 'alpha_cut_pure': np.nan}
        if r.peak_id:
            f = fits.loc[(r.cohort, r.peak_id)]
            rec.update({'n_effects': int(f.n_effects), 'n_sets_pre': int(f.n_sets_pre),
                        'n_sets_pure': int(f.n_sets_pure)})
        if not r.fine_mapped:
            status = 'not_tested' if r.fm_reason == 'not_tested' else 'not_fine_mapped'
            reason = {'not_tested': 'not tested in this cohort',
                      'no_peak': 'not fine-mapped · no peak in this cohort',
                      'outside_window': 'not fine-mapped · outside the fine-mapping window'}[r.fm_reason]
            rec.update({'cs_status': status, 'cs_reason': reason, 'cs_reason_short': reason})
            out.append(rec)
            continue
        n_eff = int(rec['n_effects'])
        if (r.cohort, r.peak_id, r.variant_id) in vkey.index:
            vv = vkey.loc[(r.cohort, r.peak_id, r.variant_id)]
            rec['alpha_pure'], rec['alpha_cut_pure'] = vv.alpha_pure, vv.alpha_cut_pure
        if r.in_cs:
            status = 'in_cs'
            reason = (f'in CS · set of {int(r.cs_size)}' if r.n_cs == 1
                      else f'in CS #{int(r.cs_index)} of {int(r.n_cs)} · set of {int(r.cs_size)}')
            short = None
        elif r.cohort_has_cs:
            status = 'below_cut'
            if r.n_cs == 1:
                reason = f'not in CS · PIP {sig4(r.pip)} < cut {sig4(r.cs_cut)}'
            else:
                cuts = ', '.join(f'#{k} {sig4(float(c))}' for k, c in (x.split(':') for x in r.cs_cuts_all.split(';')))
                reason = f'not in CS · PIP {sig4(r.pip)} · below the cut of every set ({cuts})'
            short = reason
        elif n_eff == 0:
            status, reason, short = 'no_cs_null', 'no CS · null fit (all PIP = 0)', 'no CS · null fit'
        else:
            s = top.loc[(r.cohort, r.peak_id)]
            rec['dropped_set_size'], rec['dropped_set_min_r'] = int(s.set_size), float(s.min_abs_corr)
            status, short = 'no_cs_impure', 'no CS · impure'
            dropped = f'n = {int(s.set_size)}, min |r| = {s.min_abs_corr:.2f}'
            reason = (f'no CS · impure set ({dropped})' if n_eff == 1
                      else f'no CS · {n_eff} effects fitted, no pure set (strongest: {dropped})')
        rec.update({'cs_status': status, 'cs_reason': reason, 'cs_reason_short': short or reason})
        out.append(rec)
    return pd.concat([assoc.reset_index(drop=True), pd.DataFrame(out)], axis=1)


def check_cs_status(st, eff, pairs, cs, min_abs_corr):
    """Every non-member row explained; the status partition is consistent with the fits."""
    fits = eff.drop_duplicates(['cohort', 'locus_id'])
    with_set = fits[fits.n_sets_pure >= 1]
    null = fits[fits.n_effects == 0]
    impure = fits[(fits.n_effects > 0) & (fits.n_sets_pure == 0)]
    n_no_peak = int((pairs.peak_id == '').sum())
    assert len(with_set) + len(null) + len(impure) == len(fits)
    assert len(fits) + n_no_peak == len(pairs)
    counts = st.cs_status.value_counts().to_dict()
    assert counts.get('in_cs', 0) == len(cs) and sum(counts.values()) == len(st), counts
    assert (st.in_cs == (st.cs_status == 'in_cs')).all()
    assert (st.cs_reason.fillna('') != '').all() and (st.cs_reason_short.fillna('') != '').all()
    fm = st[st.fine_mapped]
    assert (fm.pip_fitted == (fm.n_effects > 0)).all()
    assert (fm.cohort_has_cs == (fm.n_sets_pure >= 1)).all()
    assert (fm.n_cs == fm.n_sets_pure).all()
    bc = st[st.cs_status == 'below_cut']
    # the display "PIP x < cut y" depends on this; it holds on the alpha scale by
    # construction and on the PIP scale for single-effect fits
    assert (bc.alpha_pure <= bc.alpha_cut_pure).all()
    assert (bc.pip < bc.cs_cut).all(), bc[['variant_id', 'cohort', 'pip', 'cs_cut']].to_string()
    # every impure fit shows its strongest effect's dropped set, which failed the purity threshold
    si = st[st.cs_status == 'no_cs_impure'].drop_duplicates(['locus_id', 'cohort'])
    assert (si.dropped_set_min_r < min_abs_corr).all() and si.dropped_set_size.notna().all()
    # every stored set is reproduced by a surviving, pure effect (dedup: one row per distinct set)
    kept = eff[eff.stored & eff.active].drop_duplicates(['cohort', 'locus_id', 'set_size', 'min_abs_corr'])
    n_stored = int(fits.n_credible_sets_stored.sum())
    n_sets = int(cs.groupby(['cohort', 'peak_id', 'cs_index']).ngroups)
    assert len(kept) == n_stored == n_sets and (kept.purity_pass == True).all()   # noqa: E712
    dropped_single = eff[eff.active & (eff.n_effects == 1) & (eff.n_sets_pure == 0)]
    return {
        'n_fits': len(fits), 'n_pairs_without_peak': n_no_peak,
        'fits_with_set': len(with_set), 'fits_null': len(null), 'fits_impure': len(impure),
        'fits_impure_single_effect': int((impure.n_effects == 1).sum()),
        'fits_impure_multi_effect': int((impure.n_effects > 1).sum()),
        'fits_null_list': sorted(f'{r.cohort}:{r.locus_id}' for _, r in null.iterrows()),
        'dropped_single_sets_containing_lead': int(dropped_single.contains_lead.sum()),
        'n_credible_sets': n_sets, 'status_counts': counts,
        'max_abs_dpip': float(fits.max_abs_dpip.max()),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--effects', nargs='+', required=True)
    ap.add_argument('--fit-variants', nargs='+', required=True)
    ap.add_argument('--fit-json', nargs='+', required=True)
    ap.add_argument('--assoc', required=True)
    ap.add_argument('--pairs', required=True)
    ap.add_argument('--cs-members', required=True)
    ap.add_argument('--upstream-manifest', required=True)
    ap.add_argument('--tol', type=float, default=1e-3)
    ap.add_argument('--out-status', required=True)
    ap.add_argument('--out-sets', required=True)
    ap.add_argument('--out-json', required=True)
    a = ap.parse_args()
    manifest = load_json(a.upstream_manifest)
    pairs = read_tsv(a.pairs)
    cs = read_tsv(a.cs_members)
    cs['cs_index'] = num(cs.cs_index).astype(int)
    assoc = read_table(a.assoc)
    for c in ('peak_id', 'peak_tier', 'fm_reason', 'cs_cuts_all'):
        assoc[c] = assoc[c].fillna('').astype(str)
    for c in ('in_cs', 'fine_mapped', 'pip_fitted', 'cohort_has_cs'):
        assoc[c] = as_bool(assoc[c])
    eff, var, _meta = load_fits(a.effects, a.fit_variants, a.fit_json, pairs, manifest, a.tol)
    st = derive_cs_status(assoc, eff, var)
    vl = pairs.drop_duplicates('locus_id').assign(variant_id=lambda d: d.variants.str.split(',')).explode('variant_id')
    st['locus_id'] = st.variant_id.map(dict(zip(vl.variant_id, vl.locus_id)))
    assert st.locus_id.notna().all()
    summary = check_cs_status(st, eff, pairs, cs, manifest['susie']['min_abs_corr'])
    write_tsv(st, a.out_status)
    write_tsv(eff, a.out_sets)
    with open(a.out_json, 'w') as fh:
        json.dump(summary, fh, indent=2)
    log(f'CS status: {summary["status_counts"]}')


if __name__ == '__main__':
    main()
