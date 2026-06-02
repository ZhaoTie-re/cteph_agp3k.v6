import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import sys
from matplotlib.ticker import FuncFormatter

# Math symbols (italic, matching the pre-check notebook).
S_MINAC = r'$\mathit{S}_{\mathit{minAC}}$'
D_MEAN = r'$\mathit{D}_{\mathit{mean}}$'

# Colorblind-safe accents (Wong 2011).
COUNT_COLOR = '#333333'
RHO_COLOR = '#0072B2'    # blue
BETA_COLOR = '#D55E00'   # vermillion (depth effect)

INT_FMT = FuncFormatter(lambda x, _: format(int(x), ','))


def setup_style():
    plt.rcParams.update({
        'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
        'savefig.dpi': 600, 'savefig.bbox': 'tight', 'figure.dpi': 150,
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'mathtext.fontset': 'dejavusans',
        'font.size': 12, 'axes.titlesize': 14, 'axes.labelsize': 13,
        'xtick.labelsize': 11, 'ytick.labelsize': 11, 'legend.fontsize': 10.5,
        'axes.linewidth': 1.0, 'axes.edgecolor': 'black',
        'xtick.direction': 'out', 'ytick.direction': 'out',
        'axes.grid': True, 'grid.color': '#DDDDDD', 'grid.linewidth': 0.8,
        'axes.spines.top': False, 'axes.spines.right': False,
    })


def _despine(ax):
    ax.set_axisbelow(True)
    ax.grid(True, which='major')
    ax.tick_params(direction='out', length=4, width=1.0)


def plot_trends(args):
    try:
        df = pd.read_csv(args.qc_stats, sep='\t')
        required_cols = {
            'MinAC_Threshold', 'Variant_Count',
            'Spearman_Rho_SMinAC_MeanDP', 'Spearman_P_SMinAC_MeanDP',
            'Poisson_MeanDP_Beta', 'Poisson_MeanDP_95CI_Lower',
            'Poisson_MeanDP_95CI_Upper', 'Poisson_MeanDP_P',
        }
        missing_cols = required_cols - set(df.columns)
        if missing_cols:
            raise ValueError(f'qc stats file missing required columns: {sorted(missing_cols)}')

        if not args.include_zero:
            df = df[df['MinAC_Threshold'] != 0]
        df = df.sort_values('MinAC_Threshold')

        if df.empty:
            print('No data to plot (check --include-zero or input file).', file=sys.stderr)
            fig, ax = plt.subplots(figsize=(4, 4))
            ax.text(0.5, 0.5, 'No data available', ha='center', va='center')
            ax.axis('off')
            fig.savefig(args.out_pdf)
            return

        # Significance logic ---------------------------------------------------
        df['Spearman_Sig'] = df['Spearman_P_SMinAC_MeanDP'].fillna(1.0) < 0.05
        if args.sample_n:
            sample_n = args.sample_n
        elif 'Sample_Count' in df.columns:
            sample_n = int(df['Sample_Count'].iloc[0])
        else:
            sample_n = None
        if sample_n:
            threshold = 0.05 / sample_n
        else:
            threshold = 0.05
            print('Warning: no sample size; using 0.05 threshold.', file=sys.stderr)
        df['Poisson_Sig'] = df['Poisson_MeanDP_P'].fillna(1.0) < threshold

        if sample_n:
            thr_sci = f'{threshold:.1e}'
            sig_lab, nonsig_lab = rf'$p < {thr_sci}$', rf'$p \geq {thr_sci}$'
            thr_text = f'Bonferroni: 0.05 / {sample_n:,} ≈ {thr_sci}'
        else:
            sig_lab, nonsig_lab = r'$p < 0.05/N$', r'$p \geq 0.05/N$'
            thr_text = 'Significance: 0.05 / N'

        x = df['MinAC_Threshold']
        xs = sorted(df['MinAC_Threshold'].unique())

        fig, (ax_count, ax_rho, ax_depth) = plt.subplots(3, 1, figsize=(9, 11), sharex=True)

        # (a) Variant count ----------------------------------------------------
        ax_count.plot(x, df['Variant_Count'], color=COUNT_COLOR, linewidth=2, zorder=1)
        ax_count.scatter(x, df['Variant_Count'], color=COUNT_COLOR, s=42, zorder=2,
                         edgecolor='white', linewidth=0.6)
        ax_count.set_ylabel('Variant count')
        ax_count.set_title('Variant count vs MinAC threshold', loc='left', fontweight='bold', pad=10)
        ax_count.yaxis.set_major_formatter(INT_FMT)
        _despine(ax_count)

        # (b) Spearman rho (depth–burden correlation; plain p < 0.05) ----------
        ax_rho.plot(x, df['Spearman_Rho_SMinAC_MeanDP'], color=RHO_COLOR, linewidth=1.6, alpha=0.7, zorder=1)
        sig = df['Spearman_Sig']
        ax_rho.scatter(x[sig], df.loc[sig, 'Spearman_Rho_SMinAC_MeanDP'],
                       facecolors='white', edgecolors=RHO_COLOR, s=72, linewidth=1.6,
                       label=r'$p < 0.05$', zorder=2)
        ax_rho.scatter(x[~sig], df.loc[~sig, 'Spearman_Rho_SMinAC_MeanDP'],
                       color=RHO_COLOR, s=72, label=r'$p \geq 0.05$', zorder=2)
        ax_rho.axhline(0, color='#888888', linewidth=1.0)
        ax_rho.set_ylabel(rf"Spearman $\rho$  ({S_MINAC} vs {D_MEAN})")
        ax_rho.set_title('Correlation bias vs MinAC threshold', loc='left', fontweight='bold', pad=10)
        ax_rho.legend(loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0, frameon=False)
        _despine(ax_rho)

        # (c) Depth effect size (adjusted for phenotype) -----------------------
        ax_depth.fill_between(x, df['Poisson_MeanDP_95CI_Lower'], df['Poisson_MeanDP_95CI_Upper'],
                              color=BETA_COLOR, alpha=0.15, label='95% CI', zorder=0)
        ax_depth.plot(x, df['Poisson_MeanDP_Beta'], color=BETA_COLOR, linewidth=1.6, alpha=0.7, zorder=1)
        psig = df['Poisson_Sig']
        ax_depth.scatter(x[psig], df.loc[psig, 'Poisson_MeanDP_Beta'],
                         facecolors='white', edgecolors=BETA_COLOR, s=72, linewidth=1.6, label=sig_lab, zorder=2)
        ax_depth.scatter(x[~psig], df.loc[~psig, 'Poisson_MeanDP_Beta'],
                         color=BETA_COLOR, s=72, label=nonsig_lab, zorder=2)
        ax_depth.axhline(0, color='#888888', linewidth=1.0)
        ax_depth.set_ylabel(rf"Poisson $\beta_{{depth}}$  ({S_MINAC} ~ {D_MEAN})")
        ax_depth.set_title('Depth effect size vs MinAC threshold', loc='left', fontweight='bold', pad=10)
        ax_depth.text(0.02, 0.06, thr_text, transform=ax_depth.transAxes, fontsize=10,
                      bbox=dict(facecolor='white', alpha=0.95, edgecolor='#CCCCCC', boxstyle='round,pad=0.5'))
        ax_depth.legend(loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0, frameon=False)
        _despine(ax_depth)

        # Show every MinAC tick on all panels; panel letters (notebook style).
        axes = (ax_count, ax_rho, ax_depth)
        for ax in axes:
            ax.set_xticks(xs)
            ax.set_xticklabels([str(int(v)) for v in xs])
            ax.tick_params(labelbottom=True)
        ax_depth.set_xlabel('MinAC threshold')
        for ax, letter in zip(axes, 'abc'):
            ax.text(-0.06, 1.05, f'({letter})', transform=ax.transAxes,
                    fontsize=15, fontweight='bold', ha='right', va='bottom')

        fig.tight_layout()
        fig.savefig(args.out_pdf)
        print(f'Generated trend plots: {args.out_pdf}')

    except Exception as e:
        print(f'Error plotting trends: {e}', file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description='Plot MinAC trend report from merged QC summary table.')
    parser.add_argument('--qc-stats', required=True)
    parser.add_argument('--out-pdf', required=True)
    parser.add_argument('--include-zero', action='store_true', help='Include MinAC=0 in plots')
    parser.add_argument('--sample-n', type=int, help='Explicit sample size for significance correction')
    args = parser.parse_args()

    if args.sample_n is not None and args.sample_n <= 0:
        print('Error: --sample-n must be a positive integer.', file=sys.stderr)
        sys.exit(2)

    setup_style()
    plot_trends(args)


if __name__ == '__main__':
    main()
