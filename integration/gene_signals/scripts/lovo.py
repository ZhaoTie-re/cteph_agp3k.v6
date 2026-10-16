#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Leave-one-variant-out (LOVO) for the selected genes of one cohort x
#           variant set: the upstream gene test re-run on the gene's variants
#           with each variant left out in turn, so the variants that drive a
#           gene-based signal can be read off. Descriptive only: no multiple-
#           testing correction is applied to the re-runs.
#             sets     the rvtest set file: per gene one set with all its tested
#                      variants (the baseline) and one set per variant without it,
#                      written exactly as upstream writes them (one 1-bp range per
#                      variant, the contig string of the VCF), and the key table
#             collect  parse the CMC (+ cmcWald) and SKAT-O outputs and assert:
#                      every set returned once, on every sample; each set's
#                      engine count of variants equals its ranges (a left-out
#                      variant is really out); the baseline equals the published
#                      result (P, rho, beta, se, variant count) of the gene
#           The test inputs are the upstream ones: the cohort call set exported to
#           the same VCF with the same plink2 build, the published phenotype and
#           covariate files, the same rvtest options. The baseline assert is the
#           proof that the model was reproduced.
# Inputs  : see the subcommands
# Outputs : lovo.set.txt + lovo.key.tsv (sets); lovo.<cohort>.<set>.tsv (collect)
# ---------------------------------------------------------------------------
import argparse
import gzip
import json
import math
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gs_common import log, read_tsv, write_tsv  # noqa: E402

PREFIX = ('Range', 'RANGE', 'N_INFORMATIVE', 'NumVar', 'NumPolyVar')
TAIL = {'cmc': ('NonRefSite', 'Pvalue'), 'skato': ('Q', 'rho', 'Pvalue'),
        'wald': ('NonRefSite', 'Beta', 'SE', 'Pvalue')}
SEP = '__'        # set label = <set_name>__all | <set_name>__d<k>; rvtest splits names on tab / space only


def vcf_contig(path):
    """The contig string of the VCF's first record (the set file must use it verbatim)."""
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if not line.startswith('#'):
                return line.split('\t', 1)[0]
    raise SystemExit(f'ABORT: {path} has no record')


def cmd_sets(a):
    genes = read_tsv(a.genes)
    m = read_tsv(a.map)
    m = m[m.set_name.isin(genes.set_name)].copy()
    m['pos'] = m.pos.astype(int)
    bim = pd.read_csv(a.bim, sep='\t', header=None, dtype=str, usecols=[0, 1, 3], names=['chrom', 'id', 'pos'])
    n_list = sum(1 for _ in open(a.variants))
    assert len(bim) == n_list and bim.id.is_unique, (len(bim), n_list, 'the test VCF must hold every listed variant once')
    in_vcf = set(bim.id)
    first = vcf_contig(a.vcf)
    lead = 'chr' if first.startswith('chr') else ''
    key, lines = [], []
    for sn, g in m.groupby('set_name', sort=False):
        g = g.sort_values(['pos', 'variant_id'])
        assert g.variant_id.is_unique, (sn, 'a variant twice in one set')
        assert set(g.variant_id) <= in_vcf, (sn, sorted(set(g.variant_id) - in_vcf))
        assert g.chrom.nunique() == 1, sn
        # removal works through the 1-bp ranges: two variants of a gene at one position could not be separated
        assert g.pos.is_unique, (sn, 'two variants at one position; a range cannot leave one of them out')
        assert all(c not in sn for c in ' \t,:') and SEP not in sn, sn
        ranges = [f'{lead}{c}:{p}-{p}' for c, p in zip(g.chrom, g.pos)]
        vids = g.variant_id.tolist()
        n = len(vids)
        lines.append(f'{sn}{SEP}all\t' + ','.join(ranges))
        key.append({'label': f'{sn}{SEP}all', 'set_name': sn, 'gene_symbol': g.gene_symbol.iloc[0], 'k': 0,
                    'dropped': '', 'n_var': n, 'n_left': n})
        for k in range(n):
            lines.append(f'{sn}{SEP}d{k + 1}\t' + ','.join(r for i, r in enumerate(ranges) if i != k))
            key.append({'label': f'{sn}{SEP}d{k + 1}', 'set_name': sn, 'gene_symbol': g.gene_symbol.iloc[0],
                        'k': k + 1, 'dropped': vids[k], 'n_var': n, 'n_left': n - 1})
    assert lines, 'no selected gene is tested in this cohort and set'
    Path(a.out_set).write_text('\n'.join(lines) + '\n')
    write_tsv(pd.DataFrame(key), a.out_key)
    log(f'{a.cohort} {a.stratum}: {m.set_name.nunique()} genes, {len(lines)} sets '
        f'({len(lines) - m.set_name.nunique()} with one variant left out); contig "{first}"')


def read_assoc(path, method):
    """label -> dict of the rvtest columns; for cmcWald the first row (the burden term) of each block."""
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        assert tuple(head[:5]) == PREFIX and tuple(head[5:]) == TAIL[method], (path, head)
        out, order, heights = {}, [], {}
        for raw in fh:
            f = raw.rstrip('\n').split('\t')
            if f == ['']:
                continue
            assert len(f) == len(head), (path, f[:3])
            lab = f[0]
            if method == 'wald':
                if lab not in heights:
                    heights[lab] = 0
                    order.append(lab)
                    out[lab] = dict(zip(head, f))
                else:
                    assert order[-1] == lab, (path, lab, 'cmcWald blocks must be contiguous')
                heights[lab] += 1
            else:
                assert lab not in out, (path, lab, 'returned twice')
                out[lab] = dict(zip(head, f))
    return out, heights


def same(x, y):
    """Two numeric strings as rvtest / the upstream tables wrote them: equal as numbers (NA = NA)."""
    na = ('NA', 'nan', '')
    if x in na or y in na:
        return x in na and y in na
    return float(x) == float(y)


def cmd_collect(a):
    key = read_tsv(a.key)
    for c in ('k', 'n_var', 'n_left'):
        key[c] = key[c].astype(int)
    coding = json.load(open(a.coding))
    n_samples = int(coding['n_samples'])
    res = {m: read_assoc(p, m) for m, p in (('cmc', a.cmc), ('skato', a.skato), ('wald', a.wald))}
    labels = set(key.label)
    for m, (r, h) in res.items():
        assert set(r) == labels, (m, sorted(labels ^ set(r))[:5])
        if m == 'wald':
            assert set(h.values()) == {1 + a.n_covar}, (m, set(h.values()), 'cmcWald writes 1 + n_covar rows per set')
    scan = read_tsv(a.gene_scan)
    scan = scan[(scan.cohort == a.cohort) & (scan.stratum == a.stratum)].set_index(['set_name', 'method'])
    rows, n_base = [], 0
    for _, k in key.iterrows():
        for m in ('cmc', 'skato'):
            x = res[m][0][k.label]
            w = res['wald'][0][k.label] if m == 'cmc' else None
            assert int(x['N_INFORMATIVE']) == n_samples, (k.label, m, x['N_INFORMATIVE'], n_samples)
            # the engine used exactly the ranges given: a left-out variant is out
            assert int(x['NumVar']) == int(k.n_left), (k.label, m, x['NumVar'], k.n_left)
            if w is not None:
                assert (w['NumVar'], w['N_INFORMATIVE']) == (x['NumVar'], x['N_INFORMATIVE']), (k.label, 'wald')
            rec = {'cohort': a.cohort, 'stratum': a.stratum, 'set_name': k.set_name, 'gene_symbol': k.gene_symbol,
                   'method': m, 'k': int(k.k), 'dropped': k.dropped, 'n_var': int(k.n_var), 'n_left': int(k.n_left),
                   'n_var_engine': int(x['NumVar']), 'n_poly_engine': int(x['NumPolyVar']),
                   'pvalue': x['Pvalue'], 'rho': x.get('rho', 'NA'),
                   'beta': w['Beta'] if w else 'NA', 'se': w['SE'] if w else 'NA'}
            if int(k.k) == 0:
                # the baseline is the published test of the gene
                u = scan.loc[(k.set_name, m)]
                assert int(u.n_var_map) == int(k.n_var) and int(u.n_var_engine) == int(x['NumVar']), (k.label, m)
                assert int(u.n_informative) == int(x['N_INFORMATIVE']), (k.label, m)
                assert same(x['Pvalue'], u.pvalue), (k.label, m, x['Pvalue'], u.pvalue)
                if m == 'skato':
                    assert same(x['rho'], u.rho), (k.label, x['rho'], u.rho)
                else:
                    assert same(w['Beta'], u.beta) and same(w['SE'], u.se), (k.label, w['Beta'], u.beta, w['SE'], u.se)
                n_base += 1
            p = float(x['Pvalue']) if x['Pvalue'] not in ('NA', 'nan') else math.nan
            assert math.isnan(p) or 0 <= p <= 1, (k.label, m, p)
            rows.append(rec)
    out = pd.DataFrame(rows)
    assert n_base == 2 * int((key.k == 0).sum())
    out.to_csv(a.out, sep='\t', index=False)
    n_na = int(out.pvalue.isin(['NA', 'nan']).sum())
    log(f'{a.cohort} {a.stratum}: {int((key.k == 0).sum())} genes, baseline equal to the published P / rho / beta / '
        f'se for both statistics; {int((key.k > 0).sum())} leave-one-out sets x 2 statistics'
        + (f'; {n_na} P not estimable' if n_na else ''))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('sets')
    s.add_argument('--cohort', required=True)
    s.add_argument('--stratum', required=True)
    s.add_argument('--genes', required=True, help='genes.tsv of SELECT_GENES')
    s.add_argument('--map', required=True, help='the published gene map of the cohort and set (map.<set>.tsv.gz)')
    s.add_argument('--variants', required=True, help='the published variant list of the set (variants.<set>.txt)')
    s.add_argument('--bim', required=True, help='the bim of the exported test genotypes')
    s.add_argument('--vcf', required=True)
    s.add_argument('--out-set', required=True)
    s.add_argument('--out-key', required=True)
    c = sub.add_parser('collect')
    c.add_argument('--cohort', required=True)
    c.add_argument('--stratum', required=True)
    c.add_argument('--key', required=True)
    c.add_argument('--cmc', required=True, help='<out>.CMC.assoc')
    c.add_argument('--wald', required=True, help='<out>.CMCWald.assoc')
    c.add_argument('--skato', required=True, help='<out>.SkatO.assoc')
    c.add_argument('--gene-scan', required=True, help='upstream gene_scan_all.tsv')
    c.add_argument('--coding', required=True, help='upstream pheno_coding.json')
    c.add_argument('--n-covar', type=int, required=True)
    c.add_argument('--out', required=True)
    a = ap.parse_args()
    cmd_sets(a) if a.cmd == 'sets' else cmd_collect(a)


if __name__ == '__main__':
    main()
