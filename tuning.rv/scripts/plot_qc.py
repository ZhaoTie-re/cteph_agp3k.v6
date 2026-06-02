import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import argparse
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter
import gzip
import sys

# ── Publication palette (Wong 2011, colorblind-safe) + notebook aesthetic ────
GROUP_COLORS = {'Case': '#D55E00', 'PH': '#D55E00', 'Control': '#0072B2', 'NaN': '#7F7F7F'}
GROUP_FALLBACK = ['#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#F0E442']
TARGETDP_GRADIENT = ['#CFE8F3', '#73B3D8', '#2878B5', '#0F4C81']

# Math symbols (italic, matching the pre-check notebook).
S_MINAC = r'$\mathit{S}_{\mathit{minAC}}$'
D_MEAN = r'$\mathit{D}_{\mathit{mean}}$'

INT_FMT = FuncFormatter(lambda x, _: format(int(x), ','))


def setup_style():
    """White-background, despined, Arial sans-serif academic style."""
    sns.set_theme(style='ticks', context='paper')
    plt.rcParams.update({
        'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
        'savefig.dpi': 600, 'savefig.bbox': 'tight', 'figure.dpi': 150,
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'mathtext.fontset': 'dejavusans',
        'font.size': 12, 'axes.titlesize': 14, 'axes.labelsize': 13,
        'axes.titleweight': 'bold',
        'xtick.labelsize': 11, 'ytick.labelsize': 11, 'legend.fontsize': 11,
        'axes.linewidth': 1.0, 'axes.edgecolor': 'black',
        'xtick.direction': 'out', 'ytick.direction': 'out',
        'axes.grid': True, 'grid.color': '#DDDDDD', 'grid.linewidth': 0.8, 'grid.alpha': 1.0,
        'axes.spines.top': False, 'axes.spines.right': False,
    })


def _target_palette(levels):
    """Ordered (numeric-aware) sequential palette for TargetDP levels; NaN -> grey."""
    s = pd.Series(levels)
    num = s.str.extract(r'(\d+(?:\.\d+)?)', expand=False).astype(float)
    order = pd.DataFrame({'lvl': levels, 'num': num})
    order['_na'] = order['num'].isna()
    ordered = order.sort_values(['_na', 'num', 'lvl'])['lvl'].tolist()
    pal = {}
    for i, lvl in enumerate(ordered):
        pal[lvl] = '#7F7F7F' if str(lvl).lower() == 'nan' else TARGETDP_GRADIENT[i % len(TARGETDP_GRADIENT)]
    return pal, ordered


def _group_palette(levels):
    pal, fb = {}, 0
    for g in levels:
        if g in GROUP_COLORS:
            pal[g] = GROUP_COLORS[g]
        else:
            pal[g] = GROUP_FALLBACK[fb % len(GROUP_FALLBACK)]
            fb += 1
    return pal


def _despine(ax):
    ax.set_axisbelow(True)
    ax.grid(True, which='major')
    sns.despine(ax=ax)
    ax.tick_params(direction='out', length=4, width=1.0)


def load_data(args):
    try:
        df = pd.read_csv(args.sample_metrics, sep='\t')
        required_cols = {'TargetDP', 'MeanDP', 'SMinAC'}
        missing_cols = required_cols - set(df.columns)
        if missing_cols:
            raise ValueError(f'sample metrics missing required columns: {sorted(missing_cols)}')

        # Variant count = lines minus header (streamed; never loads the file).
        variant_count = 0
        with gzip.open(args.variant_metrics, 'rt') as f:
            for i, _ in enumerate(f):
                variant_count = i  # last 0-based index == line_count - 1 == data rows

        stats_df = pd.read_csv(args.qc_stats, sep='\t')
        if stats_df.empty:
            raise ValueError(f'qc stats file is empty: {args.qc_stats}')
        return df, variant_count, stats_df.iloc[0]
    except Exception as e:
        print(f'Error loading data: {e}', file=sys.stderr)
        sys.exit(1)


def get_annotation_text(stats, variant_count, sample_count, args, page_type):
    lines = [r'$\bf{Study\ design}$',
             f'Samples ($N$): {sample_count:,}',
             f'Variants ($V$): {variant_count:,}',
             f'MinAC cutoff: {args.min_ac}', '']
    if page_type == 'depth':
        try:
            lines += [r'$\bf{Distribution\ tests}$',
                      r'Mann–Whitney $U$:',
                      f"  $p = {float(stats['MWU_P_TargetDP']):.2e}$", '',
                      r'Wasserstein ($W_1$):',
                      f"  $W_1 = {float(stats['Wasserstein_Dist_TargetDP']):.4f}$"]
        except Exception:
            lines.append('Stats unavailable')
    elif page_type == 'bias':
        try:
            lines += [r'$\bf{Bias\ metrics}$',
                      r'Spearman correlation:',
                      rf"  $\rho = {float(stats['Spearman_Rho_SMinAC_MeanDP']):.3f}$",
                      rf"  $p = {float(stats['Spearman_P_SMinAC_MeanDP']):.2e}$", '',
                      'Poisson GLM:',
                      rf'  {S_MINAC} ~ Group + $Z$({D_MEAN})',
                      rf"  $\beta_{{depth}} = {float(stats['Poisson_MeanDP_Beta']):.3f}$",
                      rf"  (95% CI {float(stats['Poisson_MeanDP_95CI_Lower']):.3f}, "
                      rf"{float(stats['Poisson_MeanDP_95CI_Upper']):.3f})",
                      rf"  $p = {float(stats['Poisson_MeanDP_P']):.2e}$"]
        except Exception:
            lines.append('Stats unavailable')
    return '\n'.join(lines)


def _stats_panel(ax, text):
    ax.axis('off')
    ax.text(0.0, 0.98, text, transform=ax.transAxes, va='top', ha='left',
            fontsize=11, linespacing=1.6,
            bbox=dict(boxstyle='round,pad=0.6', facecolor='white', edgecolor='#CCCCCC', linewidth=0.8))


def _panel_letters(axes):
    for ax, letter in zip(axes, 'abc'):
        ax.text(-0.06, 1.04, f'({letter})', transform=ax.transAxes,
                fontsize=15, fontweight='bold', ha='right', va='bottom')


def plot_page1_depth(pdf, df, stats, variant_count, args):
    fig, axes = plt.subplots(1, 3, figsize=(16, 6),
                             gridspec_kw={'width_ratios': [1, 1, 0.55], 'wspace': 0.32})
    target_col = 'TargetDP'
    if target_col not in df.columns:
        return

    plot_col = 'SMinAC' if 'SMinAC' in df.columns else 'MeanDP'
    xlabel = S_MINAC if plot_col == 'SMinAC' else D_MEAN

    df = df.copy()
    df[target_col] = df[target_col].astype(str)
    levels = sorted(df[target_col].unique())
    palette, hue_order = _target_palette(levels)

    # (a) distribution
    if args.hist_stat == 'density':
        sns.kdeplot(data=df, x=plot_col, hue=target_col, hue_order=hue_order, fill=True,
                    ax=axes[0], palette=palette, alpha=0.25, linewidth=1.6, common_norm=False)
        axes[0].set_ylabel('Density')
    else:
        sns.histplot(data=df, x=plot_col, hue=target_col, hue_order=hue_order, ax=axes[0],
                     element='step', stat='count', palette=palette, alpha=0.30, linewidth=1.6)
        axes[0].set_ylabel('Sample count')
        axes[0].yaxis.set_major_formatter(INT_FMT)
    axes[0].set_title('Distribution by target depth')
    axes[0].set_xlabel(xlabel)

    # (b) ECDF
    sns.ecdfplot(data=df, x=plot_col, hue=target_col, hue_order=hue_order, ax=axes[1],
                 palette=palette, linewidth=2.4)
    axes[1].set_title('Cumulative distribution')
    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel('Cumulative proportion')

    if plot_col == 'SMinAC':
        axes[0].xaxis.set_major_formatter(INT_FMT)
        axes[1].xaxis.set_major_formatter(INT_FMT)

    for ax in (axes[0], axes[1]):
        _despine(ax)
        leg = ax.get_legend()
        if leg is not None:
            leg.set_title('Target depth')
            leg.set_frame_on(False)

    _stats_panel(axes[2], get_annotation_text(stats, variant_count, len(df), args, 'depth'))
    _panel_letters(axes)
    fig.suptitle(f'Depth distribution · MinAC ≥ {args.min_ac}', fontsize=16, y=1.02)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def plot_page2_bias(pdf, df, stats, variant_count, args):
    fig, axes = plt.subplots(1, 3, figsize=(16, 6),
                             gridspec_kw={'width_ratios': [1, 1, 0.55], 'wspace': 0.32})
    x_col, y_col = 'MeanDP', 'SMinAC'
    df = df.copy()
    df['TargetDP'] = df['TargetDP'].astype(str)

    scatter_kw = dict(alpha=0.75, s=30, edgecolor='white', linewidth=0.3)

    # (a) coloured by target depth
    t_levels = sorted(df['TargetDP'].unique())
    t_pal, t_order = _target_palette(t_levels)
    sns.scatterplot(data=df, x=x_col, y=y_col, hue='TargetDP', hue_order=t_order,
                    ax=axes[0], palette=t_pal, **scatter_kw)
    axes[0].set_title('Burden vs depth · target depth')

    # (b) coloured by phenotype group
    if 'Group' in df.columns:
        g_levels = [str(g) for g in df['Group'].fillna('NaN').unique()]
        g_pal = _group_palette(g_levels)
        sns.scatterplot(data=df.assign(Group=df['Group'].fillna('NaN').astype(str)),
                        x=x_col, y=y_col, hue='Group', palette=g_pal, ax=axes[1], **scatter_kw)
    else:
        sns.scatterplot(data=df, x=x_col, y=y_col, ax=axes[1], color='#0072B2', **scatter_kw)
    axes[1].set_title('Burden vs depth · phenotype')

    for ax, ttl in zip((axes[0], axes[1]), ('Target depth', 'Group')):
        ax.set_xlabel(D_MEAN)
        ax.set_ylabel(S_MINAC)
        ax.yaxis.set_major_formatter(INT_FMT)
        _despine(ax)
        leg = ax.get_legend()
        if leg is not None:
            leg.set_title(ttl)
            leg.set_frame_on(False)

    _stats_panel(axes[2], get_annotation_text(stats, variant_count, len(df), args, 'bias'))
    _panel_letters(axes)
    fig.suptitle(f'Depth–burden bias · MinAC ≥ {args.min_ac}', fontsize=16, y=1.02)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='Render per-MinAC QC visualization report.')
    parser.add_argument('--sample-metrics', required=True)
    parser.add_argument('--variant-metrics', required=True)
    parser.add_argument('--qc-stats', required=True)
    parser.add_argument('--min-ac', required=True)
    parser.add_argument('--out-pdf', required=True)
    parser.add_argument('--hist-stat', choices=['count', 'density'], default='density',
                        help='Statistic for the distribution panel (count or density)')
    args = parser.parse_args()

    try:
        setup_style()
        df, variant_count, stats = load_data(args)
        with PdfPages(args.out_pdf) as pdf:
            plot_page1_depth(pdf, df, stats, variant_count, args)
            plot_page2_bias(pdf, df, stats, variant_count, args)
        print(f'Generated plots: {args.out_pdf}')
    except Exception as e:
        print(f'Error generating plots: {e}', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
