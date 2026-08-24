#!/usr/bin/env python3
# ==============================================================================
# plot_gene_detail.py — one figure per hit gene.
#
# A gene-based burden test reports ONE p-value for a whole gene, which hides the
# thing that decides whether the hit is real: how many carriers there actually
# are, and whether they are spread over the gene's variants or concentrated in
# one. STBD1 is the case in point — the top gene in every stratum and method, on
# a cumulative minor-allele count of 14 (9 case, 5 control). No p-value shows
# that; this figure does.
#
# It reads what gene_detail.py already wrote, so the numbers here and in the
# per-gene tables cannot drift apart.
#
#   (a) each variant's case and control allele frequency, along the gene
#   (b) the same variants against ToMMo — the artefact check, because a variant
#       whose CONTROL frequency departs from the population reference is a
#       genotyping difference, not a disease association
#   (c) the collapsed burden: carriers and cumulative MAF, case vs control
# ==============================================================================
import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
import plot_style as S          # noqa: E402
import figure_doc               # noqa: E402

PLOT_H = 3.5
CASE_C, CTRL_C = S.ACCENT, S.DATA


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--summary', required=True, help='<GENE>.<group>.summary.txt from gene_detail.py')
    p.add_argument('--gene', required=True)
    p.add_argument('--gene-id', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--stratum', required=True)
    p.add_argument('--method', required=True)
    p.add_argument('--out-png', required=True)
    return p.parse_args()


def read_summary(path):
    """Split gene_detail's report into its key/value head and its variant table."""
    text = Path(path).read_text().splitlines()
    kv, rows, in_table, header = {}, [], False, None
    for line in text:
        if line.startswith('=== Variant Details ==='):
            in_table = True
            continue
        if in_table:
            if not line.strip():
                continue
            parts = line.split('\t')
            if header is None:
                header = parts
            elif len(parts) == len(header):
                rows.append(parts)
            continue
        m = re.match(r'^([A-Za-z0-9_\-.]+):\s*(.*)$', line)
        if m:
            kv[m.group(1)] = m.group(2).strip()
    df = pd.DataFrame(rows, columns=header) if header and rows else pd.DataFrame()
    return kv, df


def _f(v):
    try:
        x = float(v)
        return x if np.isfinite(x) else np.nan
    except (TypeError, ValueError):
        return np.nan


def freq_panel(ax, df, xcol_a, xcol_b, labels, title):
    """Per-variant frequency, case vs control (or cohort vs reference)."""
    if not len(df):
        ax.text(0.5, 0.5, 'no variant detail', transform=ax.transAxes,
                ha='center', va='center', color=S.REFERENCE)
        ax.set_xticks([]); ax.set_yticks([])
        return
    x = np.arange(len(df))
    a = df[xcol_a].map(_f).to_numpy()
    b = df[xcol_b].map(_f).to_numpy()
    w = 0.38
    ax.bar(x - w / 2, a, w, color=CASE_C, label=labels[0], zorder=3)
    ax.bar(x + w / 2, b, w, color=CTRL_C, label=labels[1], zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([s.split(':')[1] if ':' in s else s for s in df['SNPID']],
                       rotation=45, ha='right',
                       fontsize=plt.rcParams['xtick.labelsize'] - 1.5)
    ax.set_ylabel(S.EAF_SYM)
    ax.set_xlabel('variant position')
    S.panel_tag(ax, title[0])
    S.despine(ax, grid_axis='y')


def burden_panel(ax, kv):
    """The collapsed burden — the numbers the gene-level P was computed from."""
    mac_case, mac_ctrl = _f(kv.get('Cumulative_MAC_Case')), _f(kv.get('Cumulative_MAC_Ctrl'))
    maf_case, maf_ctrl = _f(kv.get('Cumulative_MAF_Case')), _f(kv.get('Cumulative_MAF_Ctrl'))
    ax.bar([0, 1], [maf_case, maf_ctrl], 0.55, color=[CASE_C, CTRL_C], zorder=3)
    for i, (v, mac) in enumerate(((maf_case, mac_case), (maf_ctrl, mac_ctrl))):
        if np.isfinite(v):
            ax.annotate(f'MAC {int(mac)}' if np.isfinite(mac) else '',
                        xy=(i, v), xytext=(0, 4), textcoords='offset points',
                        ha='center', va='bottom', fontweight='bold',
                        fontsize=plt.rcParams['legend.fontsize'] - 0.5)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['case', 'control'])
    ax.set_ylabel(S.MAF_SYM)
    top = np.nanmax([maf_case, maf_ctrl, 1e-9]) * 1.35
    ax.set_ylim(0, top)
    S.despine(ax, grid_axis='y')


def main():
    args = parse_args()
    S.setup_style()
    kv, df = read_summary(args.summary)

    fig, axes = plt.subplots(1, 3, figsize=(S.COL_DOUBLE, PLOT_H),
                             gridspec_kw={'width_ratios': [1.45, 1.45, 0.70]})
    freq_panel(axes[0], df, 'Case_AAF', 'Ctrl_AAF', ('case', 'control'),
               ('a', 'case vs control'))
    have_tommo = len(df) and 'ToMMo_AAF' in df.columns
    freq_panel(axes[1], df, 'Ctrl_AAF', 'ToMMo_AAF', ('control', 'ToMMo 60KJPN'),
               ('b', 'control vs ToMMo')) if have_tommo else \
        axes[1].axis('off')
    burden_panel(axes[2], kv)
    S.panel_tag(axes[2], 'c')
    axes[0].legend(frameon=False, fontsize=plt.rcParams['legend.fontsize'] - 1,
                   loc='upper right')

    p_method = _f(kv.get('SKAT-O_Pvalue' if args.method == 'skato' else 'CMC_Pvalue'))
    fdr = _f(kv.get('SKAT-O_FDR' if args.method == 'skato' else 'CMC_FDR'))
    n_var = kv.get('RVTest_NumVar', 'NA')
    mac_case, mac_ctrl = _f(kv.get('Cumulative_MAC_Case')), _f(kv.get('Cumulative_MAC_Ctrl'))
    orv, ci = kv.get('CMC_Effect_OR', 'NA'), kv.get('CMC_Effect_OR_95CI', '')
    stratum_h = args.stratum.replace('_', '+').upper()

    S.caption_block(
        fig, plot_h=PLOT_H,
        title=(f'{args.gene} in {args.cohort}: {n_var} variants carry a cumulative '
               f'{int(mac_case) if np.isfinite(mac_case) else "?"} case and '
               f'{int(mac_ctrl) if np.isfinite(mac_ctrl) else "?"} control minor alleles '
               f'({args.method.upper()} {S.p_tex(p_method) if np.isfinite(p_method) else "NA"}, '
               f'FDR {fdr:.3g}).' if np.isfinite(fdr) else f'{args.gene} in {args.cohort}.'),
        panels=[
            'Per-variant alternate-allele frequency in cases and controls.',
            'The same variants, control frequency against ToMMo 60KJPN — a control '
            'frequency that departs from the population reference is a genotyping '
            'difference, not an association.',
            'The collapsed burden the gene-level test acted on: cumulative minor-allele '
            'frequency over the gene, with the underlying allele counts.',
        ],
        notes=(f'{stratum_h} variants; {args.method.upper()}. '
               f'Burden OR {orv} {ci}. At this allele count the odds ratio is '
               f'estimated from a handful of carriers and its interval is wide. '
               f'Full explanation: ../README.md'),
        top_pad=0.46, wspace=0.46, left='auto', right=0.985, margin_axes=axes)
    fig.savefig(args.out_png)
    plt.close(fig)

    figure_doc.write_stats(
        args.out_png, peak=args.gene_id,
        values=[('cohort', args.cohort), ('gene', args.gene),
                ('impact stratum', args.stratum), ('test method', args.method),
                ('variants in gene', n_var),
                ('cumulative MAC, case', mac_case), ('cumulative MAC, control', mac_ctrl),
                ('cumulative MAF, case', _f(kv.get('Cumulative_MAF_Case'))),
                ('cumulative MAF, control', _f(kv.get('Cumulative_MAF_Ctrl'))),
                ('P', p_method), ('FDR', fdr),
                ('burden OR', orv), ('burden OR 95% CI', ci)])
    print(f'[plot_gene_detail] {args.cohort} {args.gene} ({args.stratum}/{args.method}): '
          f'{len(df)} variants -> {args.out_png}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in plot_gene_detail: {e}', file=sys.stderr)
        sys.exit(1)
