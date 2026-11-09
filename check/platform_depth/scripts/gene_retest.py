#!/usr/bin/env python3
"""Gene-based re-test of the Tier 1-2 genes on one genotype arm (assoc_gene's own test).

  prepare  the cohort x stratum's tested variants of the genes and their exact set-file lines
  collect  rvtest CMC / cmcWald / SKAT-O results + carriers -> one row per gene x method;
           Bonferroni from the original denominator, BH q with the gene's new P substituted into
           the original scan (every other gene unchanged); --gate checks the ORIGINAL arm against
           the published gene_scan.tsv (gate G2)

Carriers are counted at the minor allele of the ORIGINAL genotypes (one definition for every arm),
with multi-site carriers surfaced: carriers <= MAC <= 2 x n_var x carriers. No per-variant test.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bedio as B  # noqa: E402

PREFIX = ('Range', 'RANGE', 'N_INFORMATIVE', 'NumVar', 'NumPolyVar')
TAIL = {'cmc': ('NonRefSite', 'Pvalue'), 'skato': ('Q', 'rho', 'Pvalue'),
        'wald': ('NonRefSite', 'Beta', 'SE', 'Pvalue')}


def cmd_prepare(a):
    genes = a.genes.split(',')
    tested = set(Path(a.tested).read_text().split())
    lines, ids = [], []
    for g in genes:
        hit = []
        for f in sorted(Path(a.set_dir).glob('rvtest.set.chr*.txt')):
            hit += [ln for ln in f.read_text().splitlines() if ln.split('\t')[0] == g]
        if len(hit) != 1:
            sys.exit(f'ABORT: {g}: {len(hit)} set-file lines in {a.set_dir}')
        lines.append(hit[0])
        pos = {(f'chr{r.split(":")[0]}', int(r.split(':')[1].split('-')[0])) for r in hit[0].split('\t')[1].split(',')}
        ids += [v for v in tested if (v.split(':')[0], int(v.split(':')[1])) in pos]
    Path(a.out_set).write_text('\n'.join(lines) + '\n')
    Path(a.out_ids).write_text('\n'.join(sorted(set(ids))) + '\n')
    print(f'[gene_retest prepare] {len(genes)} genes, {len(set(ids))} tested variants')


def read_assoc(path, method):
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        if tuple(head[:5]) != PREFIX or tuple(head[5:]) != TAIL[method]:
            sys.exit(f'ABORT: unexpected header in {path}: {head}')
        out = {}
        for raw in fh:
            f = raw.rstrip('\n').split('\t')
            if f == ['']:
                continue
            if method == 'wald' and f[0] in out:      # cmcWald: the first row is the burden term
                continue
            if f[0] in out:
                sys.exit(f'ABORT: {f[0]} returned twice in {path}')
            out[f[0]] = dict(zip(head, f))
    return out


def bh_with(p_all, idx, p_new):
    p = p_all.copy()
    p[idx] = p_new
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return float(min(q[idx], 1.0))


def num(x):
    return float(x) if x not in ('NA', 'nan', '') else np.nan


def cmd_collect(a):
    sets = {ln.split('\t')[0]: ln.split('\t')[1] for ln in Path(a.set).read_text().splitlines() if ln.strip()}
    res = {m: read_assoc(p, m) for m, p in (('cmc', a.cmc), ('wald', a.wald), ('skato', a.skato))}
    coding = json.loads(Path(a.coding).read_text())
    pheno = pd.read_csv(a.pheno, sep='\t', dtype=str)
    cohort_ids = set(pheno.iid)                             # the minor allele is the cohort's, in every subset
    if a.include:                                          # a platform subset: only these samples count
        inc = set(Path(a.include).read_text().split())
        pheno = pheno[pheno.iid.isin(inc)]
    status = dict(zip(pheno.iid, pheno.pheno1))            # rvtest coding: 1 control, 2 case
    scan = pd.read_csv(a.gene_scan, sep='\t', dtype=str)
    den = pd.read_csv(a.denominator, sep='\t', dtype=str)
    bonf = float(den.bonferroni_threshold.iloc[0])

    # carriers at the ORIGINAL minor allele
    fam = B.read_fam(a.bfile)
    bim = B.read_bim(a.bfile)
    codes = B.unpack(B.read_rows(a.bfile, list(range(len(bim))), len(fam)), len(fam))
    ofam = B.read_fam(a.minor_from)
    obim = {f[1]: i for i, f in enumerate(B.read_bim(a.minor_from))}
    ocodes = B.unpack(B.read_rows(a.minor_from, sorted(obim.values()), len(ofam)), len(ofam))
    in_pheno = np.array([s in cohort_ids for s in ofam])
    is_case = np.array([status.get(s) == '2' for s in fam])
    is_ctrl = np.array([status.get(s) == '1' for s in fam])
    # position from the variant ID (chr<N>:<pos>:...): the .bim chromosome column has no 'chr'
    vid_pos = {f[1]: (f[1].split(':')[0], int(f[1].split(':')[1])) for f in bim}

    def minor_count(i_new, vid):
        """minor-allele count per sample in this arm (NaN = missing)."""
        oc = ocodes[sorted(obim.values()).index(obim[vid])][in_pheno]
        n5 = 2 * (oc == B.HOM5).sum() + (oc == B.HET).sum()
        n6 = 2 * (oc == B.HOM6).sum() + (oc == B.HET).sum()
        c = codes[i_new].astype(float)
        cnt = np.where(c == B.HOM5, 2, np.where(c == B.HET, 1, np.where(c == B.HOM6, 0, np.nan)))
        return cnt if n5 <= n6 else np.where(np.isnan(cnt), np.nan, 2 - cnt)

    rows = []
    for g, rng in sets.items():
        pos = {(f'chr{r.split(":")[0]}', int(r.split(':')[1].split('-')[0])) for r in rng.split(',')}
        vids = [v for v, p in vid_pos.items() if p in pos]
        idx = {f[1]: i for i, f in enumerate(bim)}
        mc = np.vstack([minor_count(idx[v], v) for v in vids])
        carried = np.nan_to_num(mc) > 0
        n_sites = carried.sum(axis=0)
        car = n_sites > 0
        multi = n_sites >= 2
        mac_case, mac_ctrl = np.nansum(mc[:, is_case]), np.nansum(mc[:, is_ctrl])
        for m in ('cmc', 'skato'):
            x = res[m].get(g)
            if x is None:
                sys.exit(f'ABORT: rvtest returned no {m} result for {g}')
            w = res['wald'].get(g) if m == 'cmc' else None
            beta, se = (num(w['Beta']), num(w['SE'])) if w else (np.nan, np.nan)
            sc = scan[(scan.cohort == a.cohort) & (scan.stratum == a.stratum) & (scan.method == m)]
            p_all = sc.pvalue.astype(float).values
            gi = np.where(sc.set_name.values == g)[0]
            if len(gi) != 1:
                sys.exit(f'ABORT: {g} is not once in the published {a.cohort} {a.stratum} {m} scan')
            p = num(x['Pvalue'])
            rows.append({
                'arm': a.arm, 'cohort': a.cohort, 'stratum': a.stratum, 'gene': g, 'method': m,
                'n_var': int(x['NumVar']), 'n_poly': int(x['NumPolyVar']), 'n_informative': int(x['N_INFORMATIVE']),
                'pvalue': p, 'rho': x.get('rho', 'NA'), 'beta': beta, 'se': se,
                'or': np.exp(beta), 'or_l95': np.exp(beta - 1.959964 * se), 'or_u95': np.exp(beta + 1.959964 * se),
                'carriers_case': int((car & is_case).sum()), 'carriers_control': int((car & is_ctrl).sum()),
                'multi_site_case': int((multi & is_case).sum()), 'multi_site_control': int((multi & is_ctrl).sum()),
                'mac_case': int(mac_case), 'mac_control': int(mac_ctrl),
                'bonferroni_threshold': bonf, 'significant_bonferroni': int(p < bonf),
                'fdr_bh': bh_with(p_all, gi[0], p), 'fdr_bh_published': float(sc.fdr_bh.values[gi[0]]),
                'pvalue_published': float(sc.pvalue.values[gi[0]])})
    out = pd.DataFrame(rows)
    out['significant_bh'] = (out.fdr_bh < 0.05).astype(int)
    for r in out.itertuples():
        nv = r.n_var
        if not (r.carriers_case + r.carriers_control <= r.mac_case + r.mac_control
                <= 2 * nv * (r.carriers_case + r.carriers_control)):
            sys.exit(f'ABORT: carrier bound violated for {r.gene} {r.method}')
    expected = len(status) if a.include else int(coding['n_samples'])
    if int(out.n_informative.iloc[0]) != expected:
        sys.exit(f'ABORT: rvtest analysed {out.n_informative.iloc[0]} samples, expected {expected}')
    out.insert(1, 'subset', a.subset)
    out.to_csv(a.out, sep='\t', index=False)
    if a.gate:
        sc = scan[(scan.cohort == a.cohort) & (scan.stratum == a.stratum)].set_index(['set_name', 'method'])
        bad = []
        for r in out.itertuples():
            pub = sc.loc[(r.gene, r.method)]
            checks = [('pvalue', r.pvalue, num(pub.pvalue)), ('n_var', r.n_var, num(pub.n_var_engine)),
                      ('carriers_case', r.carriers_case, num(pub.n_carrier_case)),
                      ('carriers_control', r.carriers_control, num(pub.n_carrier_control))]
            if r.method == 'cmc':
                checks += [('beta', r.beta, num(pub.beta)), ('se', r.se, num(pub.se))]
            for k, x, y in checks:
                if not (np.isnan(x) and np.isnan(y)) and not np.isclose(x, y, rtol=1e-5, atol=0):
                    bad.append(f'{r.gene} {r.method} {k}: {x} vs {y}')
        rep = {'cohort': a.cohort, 'stratum': a.stratum, 'rows': len(out), 'mismatches': bad, 'pass': not bad}
        Path(a.gate).write_text(json.dumps(rep, indent=2) + '\n')
        print(json.dumps(rep))
        if bad:
            sys.exit('GATE G2 FAILED')
    print(f'[gene_retest collect] {a.arm} {a.cohort} {a.stratum}: {len(out)} rows')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--genes', required=True)
    p.add_argument('--set-dir', required=True)
    p.add_argument('--tested', required=True)
    p.add_argument('--out-set', required=True)
    p.add_argument('--out-ids', required=True)
    p = sub.add_parser('collect')
    for k in ('--arm', '--cohort', '--stratum', '--set', '--cmc', '--wald', '--skato', '--coding', '--pheno',
              '--gene-scan', '--denominator', '--bfile', '--minor-from', '--out'):
        p.add_argument(k, required=True)
    p.add_argument('--gate')
    p.add_argument('--include', help='platform subset: the sample IDs that count (rvtest --peopleIncludeFile)')
    p.add_argument('--subset', default='pooled')
    a = ap.parse_args()
    {'prepare': cmd_prepare, 'collect': cmd_collect}[a.cmd](a)


if __name__ == '__main__':
    main()
