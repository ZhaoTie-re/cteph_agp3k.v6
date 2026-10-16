#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The genes of the table and their variants.
#             genes     every gene with a cross-cohort consistency tier <= MaxTier
#                       in at least one variant set, plus (IncludeBonferroni) every
#                       gene below the Bonferroni threshold in at least one cohort x
#                       variant set x statistic -- in the primary gene family (>=
#                       MinNumVar variants, as upstream; its calls and tiers are
#                       recomputed and asserted equal to upstream) or in the
#                       sensitivity family (>= SensMinNumVar). Class: Tier 1 / 2 / 3
#                       (the best primary tier), B (primary Bonferroni in one cohort
#                       only) or S (selected by the sensitivity family only).
#             variants  the union, over cohorts and variant sets, of the variants
#                       each gene was tested on (upstream variant_stats), one row
#                       per gene x variant, with the variant sets and cohorts
#                       whose test included it.
# Inputs  : upstream tiers.tsv, gene_scan_all.tsv, gene_index_all.tsv, variant_stats_all.tsv
# Outputs : genes.tsv, variants.tsv, selection.json
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gs_common import csv_list, family_calls, log, num, read_tsv, tiers_from_calls, write_tsv  # noqa: E402

CLASS_RANK = {'Tier 1': 1, 'Tier 2': 2, 'Tier 3': 3, 'B': 4, 'S': 5}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tiers', required=True)
    ap.add_argument('--gene-scan', required=True)
    ap.add_argument('--gene-index', required=True)
    ap.add_argument('--variant-stats', required=True)
    ap.add_argument('--cohorts', required=True)
    ap.add_argument('--cohort-codes', required=True)
    ap.add_argument('--strata', required=True, help='variant sets, in the order of the table')
    ap.add_argument('--max-tier', type=int, default=3)
    ap.add_argument('--include-bonferroni', default='true')
    ap.add_argument('--alpha', type=float, required=True, help='the upstream alpha')
    ap.add_argument('--min-num-var', type=int, required=True, help='the upstream minimum variants per gene (primary)')
    ap.add_argument('--sens-min-num-var', type=int, required=True, help='the sensitivity analysis\' minimum')
    ap.add_argument('--out-genes', required=True)
    ap.add_argument('--out-variants', required=True)
    ap.add_argument('--out-json', required=True)
    a = ap.parse_args()
    cohorts, strata = csv_list(a.cohorts), csv_list(a.strata)
    code = dict(zip(cohorts, csv_list(a.cohort_codes)))
    inc_bonf = str(a.include_bonferroni).lower() in ('true', '1', 'yes')
    assert 1 <= a.max_tier <= 3, a.max_tier

    tiers = read_tsv(a.tiers)
    tiers['tier'] = num(tiers.tier).astype(int)
    assert set(tiers.stratum) <= set(strata), set(tiers.stratum) - set(strata)
    scan = read_tsv(a.gene_scan)
    scan = scan[scan.cohort.isin(cohorts) & scan.stratum.isin(strata)]
    assert set(scan.cohort) == set(cohorts) and set(scan.stratum) == set(strata), 'gene scan lacks a cohort or set'
    scan['p'] = num(scan.pvalue)
    idx = read_tsv(a.gene_index)
    idx = idx[idx.cohort.isin(cohorts) & idx.stratum.isin(strata)]

    # ── the two gene families: primary (as upstream; reproduced exactly) and sensitivity ──
    prim = family_calls(scan, a.min_num_var, a.alpha)
    full = scan[num(scan.pvalue).notna()]
    assert len(prim) == len(full), 'the primary family must be every gene with a P'
    chk = prim.set_index(['cohort', 'stratum', 'method', 'set_name'])
    ref = full.set_index(['cohort', 'stratum', 'method', 'set_name'])
    assert (chk.bh == (ref.significant_bh == '1').reindex(chk.index)).all(), 'primary BH calls differ from upstream'
    assert (chk.bonf == (ref.significant_bonferroni == '1').reindex(chk.index)).all(), 'primary Bonferroni differs'
    assert (abs(chk.q - num(ref.fdr_bh).reindex(chk.index)) < 1e-6).all(), 'primary BH q differs from upstream'
    t_prim = tiers_from_calls(prim, cohorts)
    up_t = {(r.set_name, r.stratum): int(r.tier) for r in tiers.itertuples()}
    assert {k: v for k, v in t_prim.items() if v} == {k: v for k, v in up_t.items() if v}, 'primary tiers differ'
    assert a.sens_min_num_var > a.min_num_var, (a.sens_min_num_var, a.min_num_var)
    sens = family_calls(scan, a.sens_min_num_var, a.alpha)
    t_sens = tiers_from_calls(sens, cohorts)

    # ── the genes: primary tier or Bonferroni, or sensitivity tier or Bonferroni ──
    by_tier = tiers[tiers.tier <= a.max_tier]
    sel = {}
    for _, r in by_tier.iterrows():
        sel.setdefault(r.set_name, set()).add('tier')
    bonf = scan[scan.significant_bonferroni == '1']
    if inc_bonf:
        for sn in bonf.set_name.unique():
            sel.setdefault(sn, set()).add('bonferroni')
    for (sn, _st), t in t_sens.items():
        if 0 < t <= a.max_tier:
            sel.setdefault(sn, set()).add('tier_sensitivity')
    if inc_bonf:
        for sn in sens[sens.bonf].set_name.unique():
            sel.setdefault(sn, set()).add('bonferroni_sensitivity')
    assert sel, 'no gene selected'
    rows = []
    for sn, why in sel.items():
        s_idx = idx[idx.set_name == sn]
        assert len(s_idx), f'{sn} is selected but in no gene index'
        assert s_idx.gene_symbol.nunique() == 1 and s_idx.chrom.nunique() == 1, sn
        rec = {'set_name': sn, 'gene_symbol': s_idx.gene_symbol.iloc[0], 'chrom': s_idx.chrom.iloc[0],
               'pos_min': int(num(s_idx.pos_min).min()), 'pos_max': int(num(s_idx.pos_max).max()),
               'selected_by': ','.join(sorted(why))}
        best = 9
        for s in strata:
            t = tiers[(tiers.set_name == sn) & (tiers.stratum == s)]
            mapped = ((idx.set_name == sn) & (idx.stratum == s)).any()
            tv = int(t.tier.iloc[0]) if len(t) else (0 if mapped else np.nan)
            rec[f'tier_{s}'] = tv
            if tv == tv and tv > 0:
                best = min(best, tv)
            in_sens = ((sens.set_name == sn) & (sens.stratum == s)).any()
            rec[f'tier_sens_{s}'] = t_sens.get((sn, s), 0) if in_sens else np.nan
        g_scan = scan[scan.set_name == sn]
        rec['min_p'] = float(g_scan.p.min())
        rec['n_bonferroni_cells'] = int((g_scan.significant_bonferroni == '1').sum())
        rec['n_bh_cells'] = int((g_scan.significant_bh == '1').sum())
        rec['n_bonferroni_cells_sens'] = int(sens[sens.set_name == sn].bonf.sum())
        prim_sel = bool(why & {'tier', 'bonferroni'})
        rec['class'] = (f'Tier {best}' if best <= a.max_tier else 'B') if prim_sel else 'S'
        assert rec['class'] != 'B' or rec['n_bonferroni_cells'] > 0, sn
        rows.append(rec)
    genes = pd.DataFrame(rows)
    genes['_r'] = genes['class'].map(CLASS_RANK)
    genes = genes.sort_values(['_r', 'min_p', 'gene_symbol']).drop(columns='_r').reset_index(drop=True)
    genes.insert(0, 'gene_order', np.arange(1, len(genes) + 1))
    assert genes.set_name.is_unique

    # ── their variants: the union of the tested sets ────────────────────────
    vs = read_tsv(a.variant_stats)
    vs = vs[vs.set_name.isin(genes.set_name) & vs.cohort.isin(cohorts) & vs.stratum.isin(strata)]
    rows = []
    for (sn, vid), g in vs.groupby(['set_name', 'variant_id'], sort=False):
        f = g.iloc[0]
        # one annotation per variant x gene, whatever the cohort or set
        for col in ('chrom', 'pos', 'ref', 'alt', 'impact', 'effect'):
            assert g[col].nunique() == 1, (sn, vid, col, g[col].unique())
        rec = {'set_name': sn, 'variant_id': vid, 'chrom': f.chrom, 'pos': int(f.pos), 'OA': f.ref, 'EA': f.alt,
               'impact': f.impact, 'effect': f.effect}
        assert vid == f'chr{f.chrom}:{f.pos}:{f.ref}:{f.alt}', vid
        for s in strata:
            in_s = g[g.stratum == s]
            rec[f'cohorts_{s}'] = ''.join(code[c] for c in cohorts if c in set(in_s.cohort))
        rec['strata'] = ','.join(s for s in strata if rec[f'cohorts_{s}'])
        rows.append(rec)
    variants = pd.DataFrame(rows).merge(genes[['set_name', 'gene_symbol', 'gene_order']], on='set_name')
    variants = variants.sort_values(['gene_order', 'pos', 'variant_id']).reset_index(drop=True)
    # every gene tested in a cohort x set has its variants here, and the counts match the gene index
    for _, r in idx[idx.set_name.isin(genes.set_name)].iterrows():
        n = int(variants[(variants.set_name == r.set_name)][f'cohorts_{r.stratum}'].str.contains(code[r.cohort]).sum())
        assert n == int(r.n_var_map), (r.set_name, r.cohort, r.stratum, n, r.n_var_map)
    # MODERATE / HIGH variants are in both sets wherever they are tested; LOW only in the wider one
    if 'moderate_high' in strata and 'low_moderate_high' in strata:
        mh = variants[variants.impact.isin(['MODERATE', 'HIGH'])]
        assert (mh.cohorts_moderate_high == mh.cohorts_low_moderate_high).all(), 'MODERATE/HIGH set membership'
        assert (variants[variants.impact == 'LOW'].cohorts_moderate_high == '').all(), 'LOW in the narrow set'

    write_tsv(genes, a.out_genes)
    write_tsv(variants, a.out_variants)
    out = {'n_genes': len(genes), 'genes_in_order': genes.gene_symbol.tolist(),
           'class_counts': genes['class'].value_counts().to_dict(),
           'selected_by': genes.selected_by.value_counts().to_dict(),
           'n_gene_variant_rows': len(variants), 'n_variants': int(variants.variant_id.nunique()),
           'impact_counts': variants.impact.value_counts().to_dict(), 'max_tier': a.max_tier,
           'include_bonferroni': inc_bonf, 'min_num_var': a.min_num_var, 'sens_min_num_var': a.sens_min_num_var,
           'primary_reproduced': True,
           'sensitivity_tiers': {f'{sn}|{st}': t for (sn, st), t in sorted(t_sens.items()) if t}}
    json.dump(out, open(a.out_json, 'w'), indent=1)
    log(f'selection: {len(genes)} genes {out["class_counts"]}; {len(variants)} gene x variant rows '
        f'({out["n_variants"]} variants); order {", ".join(out["genes_in_order"])}')


if __name__ == '__main__':
    main()
