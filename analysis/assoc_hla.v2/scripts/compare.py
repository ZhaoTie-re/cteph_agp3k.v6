#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Put the 12 scans (3 cohorts x 2 classes x 2 models) of one case set
#           side by side.
#
#           assoc_all.tsv          every result row, stacked
#           scan_qc_all.tsv        every scan summary, stacked
#           model_concordance.tsv  per (cohort, class): fixed against random --
#                                  correlation of -log10 P and of beta over the
#                                  markers both returned, and the overlap of
#                                  their calls at the primary (GEC) threshold
#           marker_summary.tsv     one row per (class, marker): P and OR in each
#                                  of the 6 (cohort, model) cells, and for each
#                                  threshold the CROSS-COHORT CONSISTENCY of the
#                                  random (primary) model: in how many of the 3
#                                  cohorts the marker is significant with the
#                                  same direction of effect ('3/3' ... '0/3';
#                                  1/3 = "Single cohort"). `mirror_of` names the
#                                  other residue of a two-residue position whose
#                                  dosages are complements in every cohort: one
#                                  test, two rows
#           hits.tsv               the marker_summary rows significant under any
#                                  threshold in at least one cohort (random
#                                  model), or at GEC in the fixed model
#
#           The cohorts are NESTED (narrow within intermediate within full; the
#           cases are nearly the same, the controls grow), so agreement across
#           cohorts is a sample-selection sensitivity analysis, not replication.
#           The two models share samples except the relatives the fixed model
#           drops, so their agreement checks the model, not the data.
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import vocab as V                      # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assoc', nargs='+', required=True)
    p.add_argument('--qc', nargs='+', required=True)
    p.add_argument('--out-dir', required=True)
    return p.parse_args()


def cell(cohort, model):
    return f'{V.COHORT_SHORT[cohort]}_{model}'


def cohort_count(sig, orr):
    """Cohorts significant with the majority direction; (n, direction, conflict)."""
    up = int(sum(1 for s_, o in zip(sig, orr) if s_ and o > 1))
    down = int(sum(1 for s_, o in zip(sig, orr) if s_ and o < 1))
    n = max(up, down)
    direction = '' if n == 0 else ('risk' if up >= down else 'protective')
    return n, direction, bool(up and down)


def main():
    a = parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    assoc = pd.concat([pd.read_csv(f, sep='\t', dtype={'marker': str, 'chrom': str}) for f in a.assoc],
                      ignore_index=True)
    qc = pd.concat([pd.read_csv(f, sep='\t') for f in a.qc], ignore_index=True)
    order = {c: i for i, c in enumerate(V.COHORTS)}
    assoc = assoc.sort_values(['class', 'cohort', 'model', 'bp'],
                              key=lambda s: s.map(order) if s.name == 'cohort' else s).reset_index(drop=True)
    qc = qc.sort_values(['class', 'cohort', 'model'],
                        key=lambda s: s.map(order) if s.name == 'cohort' else s).reset_index(drop=True)
    assoc.to_csv(out / 'assoc_all.tsv', sep='\t', index=False, float_format='%.6g')
    qc.to_csv(out / 'scan_qc_all.tsv', sep='\t', index=False, float_format='%.6g')

    # ── fixed against random, per (cohort, class) ───────────────────────────
    rows = []
    for (cls, coh), d in assoc.groupby(['class', 'cohort'], sort=False):
        f = d[d.model == 'fixed'].set_index('marker')
        r = d[d.model == 'random'].set_index('marker')
        both = f.index.intersection(r.index)
        ok = both[f.loc[both, 'p'].notna().values & r.loc[both, 'p'].notna().values]
        lf, lr = -np.log10(f.loc[ok, 'p']), -np.log10(r.loc[ok, 'p'])
        fg, rg = set(f.index[f.sig_gec]), set(r.index[r.sig_gec])
        rows.append({'class': cls, 'cohort': coh, 'markers_both': len(ok),
                     'pearson_neglog10p': float(np.corrcoef(lf, lr)[0, 1]),
                     'spearman_p': float(stats.spearmanr(f.loc[ok, 'p'], r.loc[ok, 'p'])[0]),
                     'pearson_beta': float(np.corrcoef(f.loc[ok, 'beta'], r.loc[ok, 'beta'])[0, 1]),
                     'median_ratio_neglog10p': float(np.median(lr / lf.where(lf > 0))),
                     'gec_fixed': len(fg), 'gec_random': len(rg), 'gec_both': len(fg & rg)})
    conc = pd.DataFrame(rows)
    conc.to_csv(out / 'model_concordance.tsv', sep='\t', index=False, float_format='%.4g')

    # ── one row per (class, marker) across the 6 cells ──────────────────────
    cells = [(c, m) for c in V.COHORTS for m in V.MODELS]
    recs = []
    for (cls, mk), d in assoc.groupby(['class', 'marker'], sort=False):
        d = d.set_index(['cohort', 'model'])
        rec = {'case_set': d['case_set'].iloc[0], 'class': cls, 'marker': mk, 'gene': d['gene'].iloc[0],
               'position': d['position'].iloc[0], 'residue': d['residue'].iloc[0],
               'bp': d['bp'].iloc[0]}
        mir = d['mirror_of'].fillna('').astype(str)
        tested = [c for c in V.COHORTS if (c, V.PRIMARY_MODEL) in d.index]
        rec['mirror_of'] = mir.iloc[0] if (len(tested) == len(V.COHORTS) and mir.nunique() == 1) else ''
        for c, m in cells:
            k = cell(c, m)
            rec[f'p_{k}'] = d.loc[(c, m), 'p'] if (c, m) in d.index else np.nan
            rec[f'or_{k}'] = d.loc[(c, m), 'or'] if (c, m) in d.index else np.nan
        rec['cohorts_tested'] = len(tested)
        for model, suffix in ((V.PRIMARY_MODEL, ''), ('fixed', '_fixed')):
            for t in V.THRESHOLDS:
                if suffix and t != 'gec':
                    continue
                sig = [bool(d.loc[(c, model), f'sig_{t}']) if (c, model) in d.index else False for c in V.COHORTS]
                orr = [d.loc[(c, model), 'or'] if (c, model) in d.index else np.nan for c in V.COHORTS]
                n, direction, conflict = cohort_count(sig, orr)
                rec[f'{t}{suffix}'] = V.consistency(n)
                rec[f'{t}_n{suffix}'] = n
                if not suffix:
                    rec[f'{t}_cohorts'] = ','.join(V.COHORT_SHORT[c] for c, s_ in zip(V.COHORTS, sig) if s_) or '-'
                    if t == 'gec':
                        rec['direction'] = direction
                    rec[f'{t}_direction_conflict'] = conflict
        pr = [rec[f'p_{cell(c, V.PRIMARY_MODEL)}'] for c in V.COHORTS]
        rec['min_p'] = float(np.nanmin(pr)) if not all(np.isnan(pr)) else np.nan
        rec['max_p'] = float(np.nanmax(pr)) if not all(np.isnan(pr)) else np.nan
        ors = [rec[f'or_{cell(c, V.PRIMARY_MODEL)}'] for c in V.COHORTS]
        ors = [o for o in ors if pd.notna(o)]
        rec['same_direction'] = bool(ors) and (all(o > 1 for o in ors) or all(o < 1 for o in ors))
        recs.append(rec)
    ms = pd.DataFrame(recs)
    ms['position'] = ms['position'].astype('Int64')
    ms = ms.sort_values(['class', 'gec_n', 'min_p'], ascending=[True, False, True]).reset_index(drop=True)
    ms.to_csv(out / 'marker_summary.tsv', sep='\t', index=False, float_format='%.6g')
    anysig = (ms[[f'{t}_n' for t in V.THRESHOLDS]].sum(axis=1) > 0) | (ms['gec_n_fixed'] > 0)
    hits = ms[anysig]
    hits.to_csv(out / 'hits.tsv', sep='\t', index=False, float_format='%.6g')

    for _, r in conc.iterrows():
        print(f"[compare] {r['cohort']} {r['class']}: r(-log10P) {r['pearson_neglog10p']:.3f}, "
              f"GEC-significant fixed / random {r['gec_fixed']} / {r['gec_random']} (both {r['gec_both']})")
    for cls in V.CLASSES:
        h = ms[ms['class'] == cls]
        counts = '; '.join(f"{V.THRESHOLD_LABEL[t]}: " +
                           ', '.join(f'{n}/3 {int((h[f"{t}_n"] == n).sum())}' for n in (3, 2, 1))
                           for t in V.THRESHOLDS)
        print(f'[compare] {cls} (random model, markers): {counts}')


if __name__ == '__main__':
    main()
