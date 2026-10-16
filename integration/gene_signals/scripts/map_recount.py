#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The gene map of one cohort rebuilt BEFORE the >= MinNumVar filter, with
#           the upstream script itself (assoc_gene build_gene_map.py, staged), so
#           that the variant and gene counts the published map no longer holds
#           are exact:
#             variants  analysed; mapped / not mapped; per variant set the
#                       qualifying variants (a consequence of the set's impacts
#                       in a gene), those only in genes not tested, those in
#                       tested genes (by impact), those in > 1 tested gene
#             genes     per variant set: genes with 1, 2, >= 3 qualifying variants;
#                       variant-gene pairs in the tests
#           The rebuild must reproduce the published map exactly: the annotation
#           counts and map QC of every chromosome summed, and per variant set the
#           genes tested, genes dropped and variants tested of the denominator.
# Inputs  : the cohort's call set (PLINK) and minAC, the snpEff index, the staged
#           build_gene_map.py, the cohort's published 01.map directory
# Outputs : map_recount.<cohort>.json
# ---------------------------------------------------------------------------
import argparse
import gzip
import json
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gs_common import log, run  # noqa: E402

IMPACT_RANK = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohort', required=True)
    ap.add_argument('--bfile', required=True)
    ap.add_argument('--minac', type=int, required=True)
    ap.add_argument('--snpeff-index', required=True, help='directory of all.VQSR3.chr<N>.vcf_out.tsv.2.gz')
    ap.add_argument('--build-script', required=True, help='the staged upstream build_gene_map.py')
    ap.add_argument('--strata', nargs='+', required=True, help='<set>=<IMPACT,IMPACT> (the upstream variant sets)')
    ap.add_argument('--min-num-var', type=int, required=True)
    ap.add_argument('--map-dir', required=True, help="the cohort's published 01.map directory (cross-check)")
    ap.add_argument('--plink2', required=True)
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    c = a.cohort
    strata = {k: v.split(',') for k, v in (s.split('=', 1) for s in a.strata)}
    widest = max(strata, key=lambda s: len(strata[s]))
    assert all(set(v) <= set(strata[widest]) for v in strata.values()), 'the variant sets must be nested'
    work = Path('work_map')
    work.mkdir(exist_ok=True)

    # the analysed variants, as the upstream run listed them (plink2 --mac --write-snplist)
    run([a.plink2, '--bfile', a.bfile, '--mac', a.minac, '--write-snplist', '--threads', a.threads,
         '--out', work / 'minac'])
    keep = work / 'minac.snplist'
    by_chrom = defaultdict(list)
    with open(f'{a.bfile}.bim') as fh:
        for line in fh:
            by_chrom[line.split('\t', 1)[0]].append(line)
    chroms = sorted(by_chrom, key=lambda x: int(x) if x.isdigit() else 99)
    for ch in chroms:
        (work / f'bim.chr{ch}').write_text(''.join(by_chrom[ch]))
    del by_chrom

    def build(ch):
        idx = Path(a.snpeff_index) / f'all.VQSR3.chr{ch}.vcf_out.tsv.2.gz'
        assert idx.exists(), idx
        r = subprocess.run([sys.executable, a.build_script, '--bim', work / f'bim.chr{ch}', '--keep-variants', keep,
                            '--snpeff-index', idx, '--chrom', ch, '--impacts', *strata[widest],
                            '--gene-assignment', 'per_gene_most_severe',
                            '--out-map', work / f'map.chr{ch}.tsv.gz', '--out-anno', work / f'anno.chr{ch}.tsv',
                            '--out-qc', work / f'qc.chr{ch}.tsv'], capture_output=True, text=True)
        assert r.returncode == 0, (ch, r.stderr[-2000:])
        return ch
    with ThreadPoolExecutor(max_workers=max(1, a.threads)) as ex:
        list(ex.map(build, chroms))

    # the rebuild must be the published map
    anno = Counter()
    for ch in chroms:
        t = pd.read_csv(work / f'anno.chr{ch}.tsv', sep='\t', dtype={'Count': int}, keep_default_na=False)
        for r in t.itertuples(index=False):
            anno[(r.Impact, r.Effect, r.Biotype)] += r.Count
    pub = pd.read_csv(Path(a.map_dir) / 'anno_counts.tsv', sep='\t', dtype={'Count': int}, keep_default_na=False)
    assert anno == Counter({(r.Impact, r.Effect, r.Biotype): r.Count for r in pub.itertuples(index=False)}), \
        f'{c}: annotation counts differ from the published map'
    qc = pd.concat([pd.read_csv(work / f'qc.chr{ch}.tsv', sep='\t', dtype={'chrom': str}) for ch in chroms])
    cols = ['n_variants_in', 'n_annotated', 'n_unannotated', 'n_modifier_only', 'n_map_rows', 'n_multigene_variants']
    mine = {k: int(qc[k].sum()) for k in cols}
    tot = pd.read_csv(Path(a.map_dir) / 'map_qc.tsv', sep='\t', dtype={'chrom': str}).set_index('chrom').loc['TOTAL']
    assert mine == {k: int(tot[k]) for k in cols}, (c, mine)

    # every (variant, gene) row of the unfiltered map
    rows = []
    for ch in chroms:
        with gzip.open(work / f'map.chr{ch}.tsv.gz', 'rt') as fh:
            head = fh.readline().rstrip('\n').split('\t')
            k = {n: i for i, n in enumerate(head)}
            for line in fh:
                f = line.rstrip('\n').split('\t')
                rows.append((f[k['variant_id']], f[k['chrom']], f[k['gene_symbol']], f[k['impact']], f[k['effect']]))
    m = pd.DataFrame(rows, columns=['variant_id', 'chrom', 'gene', 'impact', 'effect'])
    del rows
    assert len(m) == mine['n_map_rows'] and m.variant_id.nunique() == mine['n_variants_in'] - mine['n_modifier_only']
    out = {'cohort': c, 'minac': a.minac, 'analysed': mine['n_variants_in'], 'mapped': int(m.variant_id.nunique()),
           'mapped_rows': len(m), 'mapped_multigene': mine['n_multigene_variants'], 'sets': {}}
    m['rank'] = m.impact.map(IMPACT_RANK)
    for s, imps in strata.items():
        x = m[m.impact.isin(imps)]
        nvar = x.groupby(['gene', 'chrom']).variant_id.nunique()
        tested = nvar[nvar >= a.min_num_var].index
        xt = x.set_index(['gene', 'chrom']).loc[tested].reset_index()
        q_best = x.sort_values(['variant_id', 'rank', 'effect']).drop_duplicates('variant_id')
        t_best = xt.sort_values(['variant_id', 'rank', 'effect']).drop_duplicates('variant_id')
        per_var = xt.groupby('variant_id').size()
        rec = {'qualifying_variants': int(x.variant_id.nunique()),
               'qualifying_by_impact': {k: int(v) for k, v in q_best.impact.value_counts().items()},
               'tested_variants': int(xt.variant_id.nunique()),
               'tested_by_impact': {k: int(v) for k, v in t_best.impact.value_counts().items()},
               'tested_multigene': int((per_var > 1).sum()),
               'tested_pairs': int(len(xt)),
               'genes': int(len(nvar)),
               'genes_1': int((nvar == 1).sum()), 'genes_2': int((nvar == 2).sum()),
               'genes_ge3': int((nvar >= 3).sum()), 'genes_tested': int(len(tested))}
        rec['untested_only_variants'] = rec['qualifying_variants'] - rec['tested_variants']
        assert all(n >= 1 for n in nvar) and rec['genes_1'] + rec['genes_2'] + rec['genes_ge3'] == rec['genes']
        # the published denominator and gene index of the set
        den = pd.read_csv(Path(a.map_dir) / f'denominator.{s}.tsv', sep='\t').iloc[0]
        gi = pd.read_csv(Path(a.map_dir) / f'gene_index.{s}.tsv', sep='\t')
        assert int(den.n_sets_chrom_split) == 0, (c, s, 'chromosome-split sets: keyed by gene and chromosome here')
        assert (rec['genes_tested'], rec['genes'] - rec['genes_tested'], rec['tested_variants'], rec['tested_pairs']) \
            == (int(den.n_genes_mapped), int(den.n_genes_dropped_lt_minvar), int(den.n_variants), int(gi.n_var_map.sum())), \
            (c, s, rec)
        assert rec['genes_2'] == int((gi.n_var_map == 2).sum()) and rec['genes_ge3'] == int((gi.n_var_map >= 3).sum())
        out['sets'][s] = rec
    json.dump(out, open(a.out, 'w'), indent=1)
    log(f'{c}: map rebuilt and equal to the published map ({mine["n_variants_in"]:,} analysed, {out["mapped"]:,} '
        f'mapped); ' + '; '.join(f'{s}: {r["qualifying_variants"]:,} qualifying variants, {r["tested_variants"]:,} in '
                                 f'{r["genes_tested"]:,} tested genes ({r["genes_2"]:,} with 2, {r["genes_ge3"]:,} with '
                                 f'>= 3), {r["untested_only_variants"]:,} only in {r["genes_1"]:,} untested genes'
                                 for s, r in out['sets'].items()))


if __name__ == '__main__':
    main()
