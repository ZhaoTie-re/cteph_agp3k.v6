#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Typing quality per PLATFORM, which is the axis the technical gradient
#           in this study runs along.
#
#           Cases and controls were sequenced differently — 3,117 controls on
#           HiSeqX at 150 bp, 327 of 452 cases on DNBSeq-T7 at 100 bp — so a
#           difference in typing quality between the groups is indistinguishable
#           from a difference in biology. It is not indistinguishable between
#           PLATFORMS, which is why every metric here is reported that way and
#           the case/control split is reported second.
#
#           Both tables are produced, by the same aggregation, and the group one is
#           a DESCRIPTION rather than a test: with zero platform overlap between the
#           groups, a difference in it has no attributable cause. It is written
#           because a reader will ask for it, and answering with a number plus that
#           warning is better than not answering.
#
#           What is measured, per platform:
#             call rate        fraction of (sample, gene) that produced two alleles
#             field resolution 2- / 3- / 4-field, the resolution actually achieved
#             ambiguity        fraction of calls where HLA-HD could not choose
#             failure          genes with no result at all — a pipeline failure
#           If these are flat across platforms, read length and depth are not
#           driving the typing. If they are not, the size of the gradient is here.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--calls', required=True)
    p.add_argument('--status', required=True)
    p.add_argument('--ambiguity', required=True)
    p.add_argument('--manifest', required=True)
    p.add_argument('--genes', required=True, help='comma-separated; the genes to score')
    p.add_argument('--out-sample', default='typing_qc_sample.tsv')
    p.add_argument('--out-platform', default='typing_qc_platform.tsv')
    p.add_argument('--out-group', default='typing_qc_group.tsv')
    p.add_argument('--out-gene', default='typing_qc_gene.tsv')
    return p.parse_args()


def depth_of(a):
    return 0 if not a or pd.isna(a) else len(str(a).split(':'))


def main():
    args = parse_args()
    calls = pd.read_csv(args.calls, sep='\t', dtype=str).fillna('')
    man = pd.read_csv(args.manifest, sep='\t', dtype=str)
    amb = pd.read_csv(args.ambiguity, sep='\t', dtype=str) if args.ambiguity else pd.DataFrame()
    genes = [g.strip() for g in args.genes.split(',') if g.strip()]

    keep = ['sample_id', 'platform', 'group', 'observed_depth']
    man = man[[c for c in keep if c in man.columns]]

    per = []
    amb_n = (amb.groupby('sample_id').size() if len(amb) else pd.Series(dtype=int))
    for _, r in calls.iterrows():
        sid = r['sample_id']
        states = [r.get(f'{g}_state', '') for g in genes]
        depths = [depth_of(r.get(f'{g}_1', '')) for g in genes
                  if r.get(f'{g}_state', '') in ('called', 'hemizygous')]
        n = len(genes)
        per.append({
            'sample_id': sid,
            'genes_scored': n,
            'called': states.count('called'),
            'hemizygous': states.count('hemizygous'),
            'not_typed': states.count('not_typed'),
            'failed': states.count('failed'),
            'call_rate': round((states.count('called') + states.count('hemizygous')) / n, 4),
            'ambiguous': int(amb_n.get(sid, 0)),
            'field2': sum(1 for d in depths if d == 2),
            'field3': sum(1 for d in depths if d == 3),
            'field4': sum(1 for d in depths if d >= 4),
            'mean_field_depth': round(sum(depths) / len(depths), 3) if depths else 0.0,
        })
    ps = pd.DataFrame(per).merge(man, on='sample_id', how='left')
    ps.to_csv(args.out_sample, sep='\t', index=False)

    num = ['call_rate', 'mean_field_depth', 'ambiguous', 'failed',
           'called', 'hemizygous', 'not_typed', 'field2', 'field3', 'field4']

    def by(col):
        """The same metric set aggregated over any one grouping column."""
        if col not in ps.columns or ps[col].isna().all():
            return pd.DataFrame(columns=[col, 'n_samples'] + [f'{m}_mean' for m in num])
        g = ps.groupby(col)[num].agg(['size', 'mean']).round(4)
        g.columns = [f'{a}_{b}' for a, b in g.columns]
        g = g.rename(columns={'call_rate_size': 'n_samples'}).reset_index()
        return g[[col, 'n_samples'] + [c for c in g.columns if c.endswith('_mean')]]

    plat = by('platform')
    plat.to_csv(args.out_platform, sep='\t', index=False)
    by('group').to_csv(args.out_group, sep='\t', index=False)

    rows = []
    for g in genes:
        col = f'{g}_state'
        if col not in calls.columns:
            continue
        vc = calls[col].value_counts()
        rows.append({'gene': g, 'called': int(vc.get('called', 0)),
                     'hemizygous': int(vc.get('hemizygous', 0)),
                     'not_typed': int(vc.get('not_typed', 0)),
                     'failed': int(vc.get('failed', 0)),
                     'call_rate': round((vc.get('called', 0) + vc.get('hemizygous', 0))
                                        / len(calls), 4)})
    pd.DataFrame(rows).to_csv(args.out_gene, sep='\t', index=False)

    print(f'[typing_qc] {len(ps)} sample(s), {len(genes)} gene(s)')
    print('    per platform (mean call rate / mean field depth / mean ambiguous / total failed):')
    for _, r in plat.iterrows():
        print(f'      {str(r["platform"]):22s} n={int(r["n_samples"]):5d}  '
              f'call {r["call_rate_mean"]:.4f}  depth {r["mean_field_depth_mean"]:.3f}  '
              f'amb {r["ambiguous_mean"]:.2f}  failed {r["failed_mean"]:.3f}')
    grp = by('group')
    if len(grp):
        print('    per group — a description, NOT a test: platform is fully confounded')
        for _, r in grp.iterrows():
            print(f'      {str(r["group"]):22s} n={int(r["n_samples"]):5d}  '
                  f'call {r["call_rate_mean"]:.4f}  depth {r["mean_field_depth_mean"]:.3f}')
    tot_fail = int(ps['failed'].sum())
    if tot_fail:
        print(f'    {tot_fail} (sample, gene) failure(s) — these are pipeline failures, '
              f'not biology, and are a differential-missingness confounder if unbalanced')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in typing_qc: {e}', file=sys.stderr)
        sys.exit(1)
