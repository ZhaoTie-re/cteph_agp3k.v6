#!/usr/bin/env python3
"""
plot_manhattan_qq.py - Publication-quality Manhattan + Q-Q plots for gene-based results.

Pipeline stage 7 (RVTEST_VISUALIZATION). Reads a post-processed association table
(containing CHR / START / Pvalue columns) and renders a two-panel figure:
  (a) Manhattan plot with a Bonferroni genome-wide threshold and labelled significant genes.
  (b) Q-Q plot with the genomic-control inflation factor (lambda_GC) and a 95% confidence band.
Writes ``<output-prefix>.png`` (600 dpi) and ``<output-prefix>.pdf`` (vector, editable text).

Usage:
    python plot_manhattan_qq.py --input <assoc> --output-prefix <prefix> [--title <str>]
"""

import argparse
import logging
import re
import sys

import matplotlib
matplotlib.use('Agg')  # Force non-interactive backend for headless (SLURM) execution
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None

logger = logging.getLogger("plot_manhattan_qq")

# ---------------------------------------------------------------------------------------
#  Publication style palette (colour-blind safe)
# ---------------------------------------------------------------------------------------
CHR_COLORS = ['#2D4F73', '#8FAFCB']   # alternating chromosome bands (dark navy / light steel)
SIG_COLOR  = '#D55E00'                # significant genes (vermillion)
THRESH_COLOR = '#C0392B'              # genome-wide threshold line (muted red)
QQ_POINT  = '#34495E'                 # Q-Q observed points (dark slate)
CI_COLOR  = '#CBD5DE'                 # 95% CI band
GRID_COLOR = '#E6E6E6'


def setup_logging() -> None:
    """Configure logging to stdout (captured into the Nextflow task log)."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Generate Manhattan and QQ plots for gene-based tests")
    parser.add_argument("--input", required=True, help="Input association file (must contain CHR, START/POS, Pvalue)")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output plot files")
    parser.add_argument("--title", required=False, default="Gene-based Test Results", help="Title for the plots")
    return parser.parse_args()


def calculate_lambda_gc(pvals) -> float:
    """Calculate the genomic-control inflation factor (lambda_GC).

    Uses the median P-value converted to a chi-squared statistic via the inverse
    survival function (precise for small P), divided by the null median chi2(1).
    """
    if len(pvals) == 0:
        return np.nan
    median_p = np.nanmedian(pvals)
    obs_median_chi2 = stats.chi2.isf(median_p, df=1)
    exp_median_chi2 = stats.chi2.ppf(0.5, df=1)
    return obs_median_chi2 / exp_median_chi2


def map_chromosome(chr_val) -> int:
    """Map a chromosome string to an integer for sorting (X=23, Y=24, XY=25, M/MT=26)."""
    s = str(chr_val).strip().lower().replace('chr', '')
    if s == 'x': return 23
    if s == 'y': return 24
    if s == 'xy': return 25
    if s == 'm' or s == 'mt': return 26
    if s.isdigit():
        return int(s)
    return 99  # Unknown/unmapped contigs -> end


def format_method_and_stratum(raw_title: str):
    """Derive a clean test-method label and (optional) impact stratum from the raw title.

    Robust to the pipeline's method tags (skato / cmc / zeggini) and severity tags such as
    'high', 'moderate_high', 'low_moderate_high' (no 'impact_' prefix required).
    """
    t = (raw_title or "").strip()
    low = t.lower()
    if "zeggini" in low:
        method = "Zeggini burden test"
    elif "cmc" in low or "burden" in low:
        method = "CMC burden test"
    elif "skato" in low or "skat-o" in low:
        method = "SKAT-O test"
    elif "skat" in low:
        method = "SKAT test"
    elif "rvtest:" in low:
        method = "Rare-variant association"
    else:
        method = t or "Rare-variant association"

    # Collect impact-severity tokens anywhere in the title (method names contain none).
    sev_map = {'low': 'LOW', 'moderate': 'MODERATE', 'mod': 'MODERATE',
               'high': 'HIGH', 'modifier': 'MODIFIER'}
    found = {sev_map[tok] for tok in re.split(r'[^a-z]+', low) if tok in sev_map}
    order = ['LOW', 'MODERATE', 'HIGH', 'MODIFIER']
    stratum = "/".join(s for s in order if s in found) if found else None
    return method, stratum


def set_pub_style() -> None:
    """Set Matplotlib rcParams for publication-quality output."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'Liberation Sans', 'sans-serif'],
        'font.size': 9,
        'axes.labelsize': 11,
        'axes.titlesize': 11,
        'xtick.labelsize': 8,
        'ytick.labelsize': 9,
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.edgecolor': '#333333',
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.size': 3.5,
        'ytick.major.size': 3.5,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'legend.fontsize': 8.5,
        'legend.frameon': False,
        'figure.dpi': 150,
        'savefig.dpi': 600,
        'pdf.fonttype': 42,   # embed editable TrueType text in the PDF/PS
        'ps.fonttype': 42,
    })


def main() -> None:
    setup_logging()
    args = parse_args()

    logger.info(f"Reading data from {args.input}...")
    try:
        df = pd.read_csv(args.input, sep='\t')
        if df.shape[1] < 2:
            df = pd.read_csv(args.input, sep=r'\s+')
    except Exception as e:
        logger.error(f"Error reading file: {e}")
        sys.exit(1)

    # Normalise column names
    df.columns = [c.lower() for c in df.columns]
    p_col = next((c for c in ['pvalue', 'p.value', 'p_value', 'p'] if c in df.columns), None)
    chr_col = next((c for c in ['chr', 'chrom'] if c in df.columns), None)
    pos_col = next((c for c in ['start', 'pos', 'position'] if c in df.columns), None)
    gene_col = next((c for c in ['gene', 'genename', 'symbol'] if c in df.columns), None)

    if not p_col or not chr_col or not pos_col:
        logger.error(f"Missing columns. Found: {df.columns.tolist()}. Need P, CHR, START.")
        sys.exit(1)

    # Keep valid rows; clamp P into (0, 1] to avoid -inf on log
    df = df.dropna(subset=[p_col, chr_col, pos_col]).copy()
    df[p_col] = df[p_col].clip(lower=1e-300, upper=1.0)
    if df.empty:
        logger.warning("No valid data rows found for plotting.")
        sys.exit(0)

    # Map chromosomes, sort, compute -log10(P) and lambda_GC
    df['CHR_NUM'] = df[chr_col].apply(map_chromosome)
    df = df[df['CHR_NUM'] < 99]
    df = df.sort_values(by=['CHR_NUM', pos_col])
    df['LOG10_P'] = -np.log10(df[p_col])
    lambda_val = calculate_lambda_gc(df[p_col].values)

    n_tests = len(df)
    bonferroni_p = 0.05 / n_tests
    bonferroni_thresh = -np.log10(bonferroni_p)
    max_logp_val = df['LOG10_P'].max()

    # =====================================================================================
    #  Figure scaffold
    # =====================================================================================
    set_pub_style()
    method, stratum = format_method_and_stratum(args.title)

    fig = plt.figure(figsize=(12.0, 4.4), facecolor='white')
    gs = fig.add_gridspec(1, 2, width_ratios=[2.6, 1.0], wspace=0.22)
    ax_man = fig.add_subplot(gs[0])
    ax_qq = fig.add_subplot(gs[1])

    # Manhattan y ceiling (leave headroom above the highest point / threshold)
    man_ylim = max(max_logp_val, bonferroni_thresh) * 1.12
    man_ylim = max(man_ylim, 6.0)

    # =====================================================================================
    #  (a) Manhattan
    # =====================================================================================
    chromosomes = sorted(df['CHR_NUM'].unique())
    ax_man.grid(axis='y', linestyle='-', linewidth=0.5, color=GRID_COLOR, zorder=0)

    # Genome-wide cumulative offsets per chromosome
    chr_offset_map = {}
    x_ticks, x_labels = [], []
    current_offset = 0
    for chrom in chromosomes:
        c_data = df[df['CHR_NUM'] == chrom]
        if c_data.empty:
            continue
        min_pos = c_data[pos_col].min()
        c_len = c_data[pos_col].max() - min_pos
        chr_offset_map[chrom] = (current_offset, min_pos)
        x_ticks.append(current_offset + c_len / 2)
        label = {23: 'X', 24: 'Y', 25: 'XY', 26: 'MT'}.get(chrom, str(chrom))
        x_labels.append(label)
        current_offset += c_len + max(c_len * 0.02, 1)  # small inter-chromosome gap

    # Background points, alternating chromosome colours (rasterised: small vector PDF)
    for i, chrom in enumerate(chromosomes):
        if chrom not in chr_offset_map:
            continue
        c_data = df[df['CHR_NUM'] == chrom]
        offset, min_p = chr_offset_map[chrom]
        x_glob = offset + (c_data[pos_col] - min_p)
        ax_man.scatter(x_glob, c_data['LOG10_P'], color=CHR_COLORS[i % 2],
                       s=8, alpha=0.9, linewidth=0, zorder=2, rasterized=True)

    # Bonferroni threshold (0.05 / number of genes)
    ax_man.axhline(bonferroni_thresh, color=THRESH_COLOR, linestyle='--', linewidth=1.0, zorder=3,
                   label=r'Bonferroni ($P = %.1e$)' % bonferroni_p)

    # Significant genes: highlight + italic labels
    sig = df[df['LOG10_P'] >= bonferroni_thresh].copy()
    if not sig.empty:
        hx, hy = [], []
        for _, row in sig.iterrows():
            c = row['CHR_NUM']
            if c not in chr_offset_map:
                continue
            offset, min_p = chr_offset_map[c]
            hx.append(offset + (row[pos_col] - min_p))
            hy.append(row['LOG10_P'])
        ax_man.scatter(hx, hy, color=SIG_COLOR, s=26, alpha=1.0, linewidth=0.4,
                       edgecolor='black', zorder=4, label='Significant')

        if gene_col:
            texts = []
            for (x, y), (_, row) in zip(zip(hx, hy), sig.iterrows()):
                texts.append(ax_man.text(x, y, str(row[gene_col]), fontstyle='italic',
                                         fontsize=8, ha='center', va='bottom', zorder=10))
            if adjust_text:
                adjust_text(texts, ax=ax_man,
                            arrowprops=dict(arrowstyle="-", color='#555555', lw=0.5),
                            expand_points=(1.4, 1.6), only_move={'text': 'xy'})

    ax_man.set_xticks(x_ticks)
    # Stagger labels onto two rows so the narrow small chromosomes (18-22) don't collide.
    ax_man.set_xticklabels([lab if i % 2 == 0 else "\n" + lab for i, lab in enumerate(x_labels)])
    ax_man.set_xlim(-current_offset * 0.01, current_offset * 1.01)
    ax_man.set_ylim(0, man_ylim)
    ax_man.set_xlabel('Chromosome', labelpad=4)
    ax_man.set_ylabel(r'$-\log_{10}(P)$', labelpad=4)
    ax_man.legend(loc='upper right', handletextpad=0.5, borderaxespad=0.4)
    ax_man.tick_params(axis='x', length=0)  # chromosome ticks read better without marks

    # =====================================================================================
    #  (b) Q-Q
    # =====================================================================================
    p_sorted = np.sort(df[p_col].values)
    observed_logp = -np.log10(p_sorted)
    n_points = len(df)
    idx = np.arange(1, n_points + 1)
    expected_logp = -np.log10((idx - 0.5) / n_points)

    # 95% confidence band (Beta order statistics)
    lower_log = -np.log10(stats.beta.ppf(0.975, idx, n_points - idx + 1))
    upper_log = -np.log10(stats.beta.ppf(0.025, idx, n_points - idx + 1))
    ax_qq.fill_between(expected_logp, lower_log, upper_log, color=CI_COLOR, alpha=0.7,
                       linewidth=0, zorder=1, label='95% CI')

    qq_xmax = expected_logp.max() * 1.05
    qq_ymax = max(observed_logp.max(), bonferroni_thresh) * 1.08
    diag = max(qq_xmax, qq_ymax)
    ax_qq.plot([0, diag], [0, diag], color=THRESH_COLOR, linestyle='--', linewidth=1.0,
               zorder=2, label='Null')
    ax_qq.scatter(expected_logp, observed_logp, c=QQ_POINT, s=12, alpha=0.85,
                  linewidth=0, zorder=3, rasterized=True)

    ax_qq.set_xlim(0, qq_xmax)
    ax_qq.set_ylim(0, qq_ymax)
    ax_qq.set_box_aspect(1)  # square panel without distorting the data ranges
    ax_qq.grid(True, linestyle='-', linewidth=0.5, color=GRID_COLOR, zorder=0)
    ax_qq.set_xlabel(r'Expected $-\log_{10}(P)$', labelpad=4)
    ax_qq.set_ylabel(r'Observed $-\log_{10}(P)$', labelpad=4)

    stats_text = f"$\\lambda_{{GC}} = {lambda_val:.3f}$\n$N = {n_points:,}$ genes"
    ax_qq.text(0.04, 0.96, stats_text, transform=ax_qq.transAxes, fontsize=9,
               va='top', ha='left',
               bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#999999', linewidth=0.8))

    # =====================================================================================
    #  Titles / panel labels
    # =====================================================================================
    suptitle = method + (f"   ·   impact: {stratum}" if stratum else "")
    fig.suptitle(suptitle, fontsize=12, fontweight='bold', y=1.02)
    for ax, lab in ((ax_man, 'a'), (ax_qq, 'b')):
        ax.text(-0.02, 1.06, lab, transform=ax.transAxes, fontsize=13, fontweight='bold',
                va='bottom', ha='right')

    fig.savefig(f"{args.output_prefix}.png", bbox_inches='tight')
    fig.savefig(f"{args.output_prefix}.pdf", bbox_inches='tight')
    logger.info(f"Saved plots to {args.output_prefix}.png/pdf  (lambda_GC={lambda_val:.3f}, N={n_points})")


if __name__ == "__main__":
    main()
