#!/usr/bin/env python3
# ==============================================================================
# cross_cohort_genes.py — every gene that reached a tier in ANY cohort, reported
# in ALL of them.
#
# The gene-level counterpart of _shared/scripts/cross_cohort.py, and it makes the
# same distinction that one does: a blank cell is never inferred. A gene can be
# absent from a cohort's table for two different reasons — it was tested and did
# not reach a tier (`not_a_hit`), or it was never tested there at all because too
# few qualifying variants survived that cohort's minAC (`not_tested`). Collapsing
# those two into one "not significant" throws away the distinction a reader most
# wants.
#
# The three cohorts are NESTED (narrow subset of intermediate subset of full) and
# share their cases entirely, so agreement between them is NOT replication. It is
# a sensitivity check on the ancestry filter. Every output here is labelled that
# way so the tables cannot be read as independent support.
# ==============================================================================
import argparse
import sys

import numpy as np
import pandas as pd

TIER_RANK = {'significant': 0, 'suggestive': 1, 'not_a_hit': 2, 'not_tested': 3}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scan', action='append', required=True,
                   metavar='COHORT=STRATUM=METHOD=PATH',
                   help='gene_scan.tsv for one (cohort, stratum, method). Repeatable.')
    p.add_argument('--qc', action='append', default=[], metavar='COHORT=STRATUM=METHOD=PATH')
    p.add_argument('--cohort-order', required=True, help='comma-separated, narrowest first')
    p.add_argument('--out-genes', default='gene_crosscohort.tsv')
    p.add_argument('--out-scans', default='scan_qc_all.tsv')
    p.add_argument('--out-hits', default='gene_hits_all.tsv')
    return p.parse_args()


def spec(entries):
    out = {}
    for e in entries:
        parts = e.split('=', 3)
        if len(parts) != 4:
            raise SystemExit(f'--scan/--qc must be COHORT=STRATUM=METHOD=PATH, got {e!r}')
        out[(parts[0], parts[1], parts[2])] = parts[3]
    return out


def main():
    args = parse_args()
    order = [c.strip() for c in args.cohort_order.split(',') if c.strip()]
    scans = spec(args.scan)

    frames = []
    for (cohort, stratum, method), path in sorted(scans.items()):
        d = pd.read_csv(path, sep='\t')
        d['cohort'], d['stratum'], d['method'] = cohort, stratum, method
        frames.append(d)
    if not frames:
        raise SystemExit('no scan tables given')
    allscan = pd.concat(frames, ignore_index=True)
    allscan.to_csv(args.out_hits, sep='\t', index=False)

    qc = spec(args.qc)
    if qc:
        qframes = []
        for (cohort, stratum, method), path in sorted(qc.items()):
            qframes.append(pd.read_csv(path, sep='\t'))
        pd.concat(qframes, ignore_index=True).to_csv(args.out_scans, sep='\t', index=False)

    rows = []
    # One comparison per (stratum, method): a gene's tier is only comparable to
    # another cohort's under the SAME test on the SAME variant stratum.
    for (stratum, method), grp in allscan.groupby(['stratum', 'method']):
        tested = {c: set(grp[grp.cohort == c]['Gene']) for c in order}
        hits = grp[grp['tier'].isin(['significant', 'suggestive'])]
        union = sorted(set(hits['Gene']))
        for gene in union:
            best = min((TIER_RANK.get(t, 9) for t in hits[hits.Gene == gene]['tier']),
                       default=9)
            best_tier = next(k for k, v in TIER_RANK.items() if v == best)
            for c in order:
                sel = grp[(grp.cohort == c) & (grp.Gene == gene)]
                if not len(sel):
                    # Never tested here — not the same as tested and null.
                    rows.append(dict(stratum=stratum, method=method, Gene=gene,
                                     best_tier=best_tier, cohort=c,
                                     called='not_tested', Pvalue=np.nan, FDR_BH=np.nan,
                                     NumVar=np.nan, chrom=np.nan))
                    continue
                r = sel.iloc[0]
                rows.append(dict(stratum=stratum, method=method, Gene=gene,
                                 best_tier=best_tier, cohort=c,
                                 called=r['tier'], Pvalue=r['Pvalue'],
                                 FDR_BH=r['FDR_BH'], NumVar=r.get('NumVar', np.nan),
                                 chrom=r.get('chrom', np.nan)))
    out = pd.DataFrame(rows)
    if len(out):
        out['_c'] = out['cohort'].map({c: i for i, c in enumerate(order)})
        out['_t'] = out['best_tier'].map(TIER_RANK)
        out = out.sort_values(['stratum', 'method', '_t', 'Gene', '_c']).drop(columns=['_c', '_t'])
    out.to_csv(args.out_genes, sep='\t', index=False)

    n_gene = out['Gene'].nunique() if len(out) else 0
    n_sig = out[out.best_tier == 'significant']['Gene'].nunique() if len(out) else 0
    print(f'[cross_cohort_genes] {len(scans)} scans, {n_gene} distinct gene(s) reaching a '
          f'tier in >=1 cohort ({n_sig} significant somewhere) -> {args.out_genes}')
    if len(out):
        shared = (out[out.called.isin(['significant', 'suggestive'])]
                  .groupby(['stratum', 'method', 'Gene']).size().value_counts().sort_index())
        for k, v in shared.items():
            print(f'    a hit in {k} of {len(order)} cohorts: {v} gene-scan pair(s)')
        print('    NOTE: the cohorts are nested and share their cases, so this counts '
              'ROBUSTNESS to the ancestry filter, not replication.')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in cross_cohort_genes: {e}', file=sys.stderr)
        sys.exit(1)
