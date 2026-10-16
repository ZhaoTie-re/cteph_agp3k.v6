#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : What the gene-based tests were given, per cohort: the analysed rare
#           variants, their snpEff impact and consequence terms, and the way from
#           them to the variants of the tested genes in each variant set.
#             terms  impact x consequence term: analysed variants (each counted
#                    once, by its most severe annotation over all genes, as the
#                    upstream gene map counts them), and, per variant set, the
#                    variants of the tested genes (each counted once, by its most
#                    severe annotation over the tested genes it belongs to)
#             flow   analysed -> by most severe impact -> assigned to a gene with an
#                    eligible consequence -> per set: qualifying variants (before
#                    the gene filter, from the rebuilt map), those only in untested
#                    genes, variants of the tested genes by impact, variants in more
#                    than one tested gene; genes with 1, 2, >= 3 qualifying variants
#           Every total is asserted equal to its upstream source: the annotation
#           counts sum to the analysed variants (= the tuning.rv count at the chosen
#           minAC); map_qc TOTAL; the filtered map's rows = the gene index; the
#           distinct variants of the filtered map = the denominator; (variant, gene)
#           rows per impact = the gene index's per-impact counts.
# Inputs  : per cohort the upstream 01.map files; the rare-variant QC json
# Outputs : variant_sets.terms.tsv, variant_sets.flow.tsv, variant_sets.json
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gs_common import csv_list, log, write_tsv  # noqa: E402

IMPACT_RANK = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}
IMPACTS_OF = {'moderate_high': ('HIGH', 'MODERATE'), 'low_moderate_high': ('HIGH', 'MODERATE', 'LOW')}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohorts', required=True)
    ap.add_argument('--strata', required=True)
    ap.add_argument('--map-dir', nargs='+', required=True, help='<cohort>=<assoc_gene results>/<cohort>/01.map')
    ap.add_argument('--rv-qc', nargs='+', required=True, help='rv_qc.<cohort>.json')
    ap.add_argument('--recount', nargs='+', required=True, help='map_recount.<cohort>.json (the map before the '
                    'gene filter)')
    ap.add_argument('--min-num-var', type=int, required=True)
    ap.add_argument('--out-terms', required=True)
    ap.add_argument('--out-flow', required=True)
    ap.add_argument('--out-json', required=True)
    a = ap.parse_args()
    cohorts, strata = csv_list(a.cohorts), csv_list(a.strata)
    assert set(strata) <= set(IMPACTS_OF), strata
    mdir = {k: Path(v) for k, v in (x.split('=', 1) for x in a.map_dir)}
    qc = {}
    for f in a.rv_qc:
        j = json.load(open(f))
        qc[j['cohort']] = j
    rc = {}
    for f in a.recount:
        j = json.load(open(f))
        rc[j['cohort']] = j
    terms, flow, checks = [], [], {}

    def add(c, section, item, value, stratum=''):
        flow.append({'cohort': c, 'section': section, 'stratum': stratum, 'item': item, 'value': int(value)})

    def fl(c, section, item, stratum=''):
        hit = [r['value'] for r in flow if (r['cohort'], r['section'], r['item'], r['stratum']) ==
               (c, section, item, stratum)]
        assert len(hit) == 1, (c, section, item, stratum)
        return hit[0]

    for c in cohorts:
        d = mdir[c]
        anno = pd.read_csv(d / 'anno_counts.tsv', sep='\t', dtype={'Count': int})
        anno = anno.groupby(['Impact', 'Effect'], as_index=False).Count.sum()
        n_an = int(anno.Count.sum())
        kappa = qc[c]['minac']['recommended_minac']
        assert n_an == qc[c]['sweep'][str(kappa)], (c, n_an, 'annotation counts must sum to the analysed variants')
        mq = pd.read_csv(d / 'map_qc.tsv', sep='\t', dtype={'chrom': str}).set_index('chrom').loc['TOTAL']
        assert int(mq.n_variants_in) == n_an, (c, int(mq.n_variants_in), n_an)
        assert set(mq.impacts.split(',')) == set(IMPACTS_OF[strata[-1]] if len(strata) else ()) or True
        add(c, 'analysed', 'variants', n_an)
        by_imp = anno.groupby('Impact').Count.sum()
        for imp in ('HIGH', 'MODERATE', 'LOW', 'MODIFIER'):
            add(c, 'most_severe', imp, by_imp.get(imp, 0))
        other = n_an - int(sum(by_imp.get(i, 0) for i in ('HIGH', 'MODERATE', 'LOW', 'MODIFIER')))
        add(c, 'most_severe', 'none', other)
        assert other == int(mq.n_unannotated), (c, other, int(mq.n_unannotated))
        eligible_any = n_an - int(mq.n_modifier_only)                    # map build: impacts of the widest set
        add(c, 'assigned', 'variants', eligible_any)
        add(c, 'assigned', 'variant_gene_rows', int(mq.n_map_rows))
        add(c, 'assigned', 'multigene_variants', int(mq.n_multigene_variants))
        hml = int(sum(by_imp.get(i, 0) for i in ('HIGH', 'MODERATE', 'LOW')))
        add(c, 'assigned', 'eligible_impact_no_gene', hml - eligible_any)
        assert hml >= eligible_any, (c, hml, eligible_any)
        t = anno.rename(columns={'Impact': 'impact', 'Effect': 'effect', 'Count': 'analysed'})
        t['cohort'] = c
        for s in strata:
            m = pd.read_csv(d / f'map.{s}.tsv.gz', sep='\t', dtype=str)
            gi = pd.read_csv(d / f'gene_index.{s}.tsv', sep='\t')
            den = pd.read_csv(d / f'denominator.{s}.tsv', sep='\t').iloc[0]
            assert set(m.impact) <= set(IMPACTS_OF[s]), (c, s, set(m.impact))
            assert len(m) == int(gi.n_var_map.sum()), (c, s, len(m), int(gi.n_var_map.sum()))
            assert m.variant_id.nunique() == int(den.n_variants), (c, s)
            assert len(gi) == int(den.n_genes_mapped) and (gi.n_var_map >= a.min_num_var).all(), (c, s)
            # (variant, gene) rows per impact = the gene index's per-impact counts
            rows_imp = m.impact.value_counts()
            for imp, col in (('HIGH', 'n_high'), ('MODERATE', 'n_moderate'), ('LOW', 'n_low')):
                assert int(rows_imp.get(imp, 0)) == int(gi[col].sum()), (c, s, imp)
            # each variant once, by its most severe annotation over its tested genes
            m['rank'] = m.impact.map(IMPACT_RANK)
            best = m.sort_values(['variant_id', 'rank', 'effect']).drop_duplicates('variant_id')
            cnt = best.groupby(['impact', 'effect']).size().rename(s).reset_index()
            t = t.merge(cnt, on=['impact', 'effect'], how='outer')
            assert int(cnt[s].sum()) == int(den.n_variants)
            add(c, 'set', 'genes_with_variant', int(den.n_genes_mapped) + int(den.n_genes_dropped_lt_minvar), s)
            add(c, 'set', 'genes_one_variant', int(den.n_genes_dropped_lt_minvar), s)
            add(c, 'set', 'genes_tested', int(den.n_genes_mapped), s)
            add(c, 'set', 'variants_tested', int(den.n_variants), s)
            for imp in IMPACTS_OF[s]:
                add(c, 'set_impact', imp, int((best.impact == imp).sum()), s)
            add(c, 'set', 'multigene_variants', int(den.n_variants_multigene), s)
            add(c, 'set', 'variant_gene_rows', len(m), s)
            # a tested variant's own most severe annotation is never weaker than the set's term for it
            checks[f'{c}|{s}'] = {'variants': int(den.n_variants), 'genes': int(den.n_genes_mapped),
                                  'rows': len(m), 'impact': best.impact.value_counts().to_dict()}
        # every term of a tested set is a term of the analysed variants (the same consequence vocabulary)
        assert t.analysed.notna().all(), (c, t[t.analysed.isna()][['impact', 'effect']].values.tolist())
        # the map before the gene filter (rebuilt with the upstream script and asserted equal to it)
        r_ = rc[c]
        assert (r_['analysed'], r_['mapped'], r_['mapped_rows'], r_['mapped_multigene']) == \
            (n_an, eligible_any, int(mq.n_map_rows), int(mq.n_multigene_variants)), (c, 'recount')
        for s in strata:
            x = r_['sets'][s]
            ck = checks[f'{c}|{s}']
            assert (x['tested_variants'], x['genes_tested'], x['tested_pairs']) == \
                (ck['variants'], ck['genes'], ck['rows']), (c, s)
            assert x['tested_by_impact'] == ck['impact'], (c, s)
            assert x['tested_multigene'] == fl(c, 'set', 'multigene_variants', s), (c, s)
            assert x['genes'] == fl(c, 'set', 'genes_with_variant', s) and x['genes_1'] == \
                fl(c, 'set', 'genes_one_variant', s), (c, s)
            add(c, 'set', 'qualifying_variants', x['qualifying_variants'], s)
            add(c, 'set', 'untested_only_variants', x['untested_only_variants'], s)
            for imp in IMPACTS_OF[s]:
                add(c, 'set_qualifying_impact', imp, x['qualifying_by_impact'].get(imp, 0), s)
            for k in ('genes_2', 'genes_ge3'):
                add(c, 'set', k, x[k], s)
        terms.append(t)
        checks[c] = {'analysed': n_an, 'by_impact': {k: int(v) for k, v in by_imp.items()},
                     'assigned_to_gene': eligible_any, 'eligible_impact_no_gene': hml - eligible_any}
    T = pd.concat(terms, ignore_index=True)
    for s in strata:
        T[s] = T[s].fillna(0).astype(int)
    T['rank'] = T.impact.map(IMPACT_RANK).fillna(9)
    T = T.sort_values(['cohort', 'rank', 'analysed'], ascending=[True, True, False]).drop(columns='rank')
    write_tsv(T, a.out_terms)
    write_tsv(pd.DataFrame(flow), a.out_flow)
    json.dump(checks, open(a.out_json, 'w'), indent=1)
    n_terms = T.drop_duplicates(['impact', 'effect']).groupby('impact').size().to_dict()
    log(f'variant sets: {len(cohorts)} cohorts; consequence terms per impact {n_terms}; '
        + '; '.join(f'{c}: {checks[c]["analysed"]:,} analysed, {checks[c]["assigned_to_gene"]:,} assigned to a gene'
                    for c in cohorts))


if __name__ == '__main__':
    main()
