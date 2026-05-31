#!/usr/bin/env python3
"""
plot_regional_comparison.py

Publication-style stacked regional association plot comparing GWAS results
across QC stages: pre-genotype QC, post-genotype QC, and post-variant QC.

Main design choices
-------------------
- All plotted variants are circular points.
- The target variant is highlighted by a red circular ring plus inner dot.
- The target variant is annotated using the value passed to --target-loci by default.
- The x-axis is compact and the y-direction is visually taller.
- Panel labels do not include a/b/c prefixes.
- Panels are separated with moderate white space for a cleaner journal-style layout.
- The style is tuned for a clean academic / Nature-like figure.

Usage
-----
    python plot_regional_comparison.py \
        --glm-files *.PHENO1.glm.logistic \
        --ld-files  ld_*.vcor \
        --target-loci "chr16:53887925:T:C" \
        --genomic-locus "chr16:53703963-54121941" \
        --gene-gtf /path/to/gencode.gtf.gz \
        --out regional_comparison \
        --pdf
"""

import argparse
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import FuncFormatter, MaxNLocator
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator


# ── Constants ─────────────────────────────────────────────────────────────────

STAGE_ORDER = ["pre_genotype_qc", "post_genotype_qc", "post_variant_qc"]

STAGE_LABELS = {
    "pre_genotype_qc": "Pre-genotype QC",
    "post_genotype_qc": "Post-genotype QC",
    "post_variant_qc": "Post-variant QC",
}

PANEL_FACE = "#FFFFFF"
FIGURE_FACE = "#FFFFFF"

NO_LD_COLOR = "#B8C2CF"
TARGET_COLOR = "#B2182B"
GENE_LINE_COLOR = "#748294"
GENE_EXON_COLOR = "#2F6EA6"

GRID_COLOR = "#DDE4EE"
SPINE_COLOR = "#1F2937"
TEXT_COLOR = "#111827"
MUTED_TEXT_COLOR = "#5B6472"

THRESHOLD_COLOR = "#9B5F5F"
TARGET_GUIDE_COLOR = "#D77272"

POINT_SIZE = 9.0

GENE_FILTER_RE = re.compile(
    r"^(?:"
    r"AC\d|AL\d|AP\d|BX\d|CT[A-Z]?\d|CU\d|FP\d|LOC\d|LINC\d|"
    r"MIR\d|SNOR[A-Z]?\d|RNU\d|RNA\d|RF\d|RP\d|RN7SL|RNF\dP|"
    r"IG[HKL][A-Z0-9-]*|TR[ABDG][A-Z0-9-]*P?"
    r")",
    re.IGNORECASE,
)

LD_BIN_EDGES = [0.0, 0.2, 0.4, 0.6, 0.8, 1.000001]
LD_BIN_LABELS = ["0.0 - 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1.0"]
LD_BIN_COLORS = [
    "#3D5EA8",
    "#3E95BE",
    "#69B494",
    "#D0AB45",
    "#C53A2D",
]

LD_CMAP = mcolors.ListedColormap(LD_BIN_COLORS)
LD_NORM = mcolors.BoundaryNorm(LD_BIN_EDGES, LD_CMAP.N)


# ── Argument parsing ──────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    p.add_argument(
        "--glm-files",
        nargs="+",
        required=True,
        help="plink2 .glm.* result files, one or more QC stages",
    )
    p.add_argument(
        "--ld-files",
        nargs="+",
        required=True,
        help="plink2 LD files, e.g. ld_<stage>.vcor/.vcor1/.ld",
    )
    p.add_argument(
        "--target-loci",
        required=True,
        help="Target variant ID, e.g. chr16:53887925:T:C",
    )
    p.add_argument(
        "--target-label",
        default=None,
        help="Optional label for the target variant. Default: use --target-loci value.",
    )
    p.add_argument(
        "--hide-target-label",
        action="store_true",
        help="Hide the text label next to the target variant.",
    )
    p.add_argument(
        "--genomic-locus",
        default=None,
        help="Optional locus string chrN:START-END; used for gene track range",
    )
    p.add_argument(
        "--gene-gtf",
        default=None,
        help="Optional hg38 GTF(.gz) path; if provided, add gene structure panel",
    )
    p.add_argument(
        "--max-track-genes",
        type=int,
        default=30,
        help="Max number of genes to draw in gene track",
    )
    p.add_argument(
        "--show-uncharacterized-genes",
        action="store_true",
        help="Show low-information gene symbols such as AC*, AL*, LINC*, MIR*, etc.",
    )
    p.add_argument(
        "--pdf",
        action="store_true",
        help="Also save vector PDF",
    )
    p.add_argument(
        "--out",
        required=True,
        help="Output filename prefix",
    )
    p.add_argument(
        "--title",
        default="Regional association QC-stage comparison",
        help="Figure title",
    )
    p.add_argument(
        "--y-max",
        type=float,
        default=None,
        help="Optional fixed y-axis maximum for -log10(P)",
    )
    p.add_argument(
        "--fig-width",
        type=float,
        default=6.25,
        help="Figure width in inches",
    )
    p.add_argument(
        "--fig-height",
        type=float,
        default=None,
        help="Optional figure height in inches",
    )

    return p.parse_args()


# ── Utility ───────────────────────────────────────────────────────────────────

def _chr_variants(variant_id: str) -> set:
    variants = {variant_id}

    if variant_id.startswith("chr"):
        variants.add(variant_id[3:])
    else:
        variants.add("chr" + variant_id)

    return variants


def _chr_aliases(chrom: str) -> set:
    c = str(chrom).strip()

    if not c:
        return set()

    if c.startswith("chr"):
        return {c, c[3:]}

    return {c, "chr" + c}


def infer_stage(filename: str):
    m = re.search(r"\.region_check\.([a-z_]+)\.", filename)
    return m.group(1) if m else None


def infer_stage_from_ld(filename: str):
    m = re.search(r"ld_([a-z_]+)\.(vcor|vcor1|ld)", os.path.basename(filename))
    return m.group(1) if m else None


def parse_genomic_locus(locus: str):
    m = re.match(r"^(chr?[0-9XYM]+):(\d+)-(\d+)$", str(locus).strip())

    if not m:
        return None

    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))

    if start > end:
        start, end = end, start

    return chrom, start, end


def parse_gtf_attr(attr: str) -> dict:
    parsed = {}

    for chunk in str(attr).strip().split(";"):
        chunk = chunk.strip()

        if not chunk or " " not in chunk:
            continue

        k, v = chunk.split(" ", 1)
        parsed[k.strip()] = v.strip().strip('"')

    return parsed


def use_publication_style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": FIGURE_FACE,
            "axes.facecolor": PANEL_FACE,
            "savefig.facecolor": FIGURE_FACE,
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
            "font.size": 8.2,
            "axes.labelsize": 9.2,
            "axes.titlesize": 9.4,
            "axes.titleweight": "semibold",
            "xtick.labelsize": 8.1,
            "ytick.labelsize": 8.1,
            "legend.fontsize": 7.5,
            "axes.linewidth": 0.78,
            "xtick.major.width": 0.70,
            "ytick.major.width": 0.70,
            "xtick.major.size": 3.2,
            "ytick.major.size": 3.2,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "axes.unicode_minus": False,
            "mathtext.default": "regular",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def make_ld_legend_handles() -> list:
    handles = []

    for color, label in zip(reversed(LD_BIN_COLORS), reversed(LD_BIN_LABELS)):
        handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="None",
                markerfacecolor=color,
                markeredgecolor="white",
                markeredgewidth=0.25,
                markersize=4.5,
                label=label,
            )
        )

    return handles


def should_plot_gene(gene_name: str, show_uncharacterized: bool = False) -> bool:
    clean_name = str(gene_name or "").strip()

    if not clean_name:
        return False

    if show_uncharacterized:
        return True

    return GENE_FILTER_RE.match(clean_name) is None


def format_mb_ticks(x, _):
    return f"{x:.2f}".rstrip("0").rstrip(".")


def format_p_value(p: float) -> str:
    if not np.isfinite(p) or p <= 0:
        return "NA"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.4f}".rstrip("0").rstrip(".")


def draw_gene_direction_arrows(ax, x0: float, x1: float, y: float, strand: str) -> None:
    span = x1 - x0

    if span <= 0.018:
        return

    n_arrows = int(np.clip(np.floor(span / 0.065), 1, 6))
    positions = np.linspace(x0 + span * 0.18, x1 - span * 0.18, n_arrows)
    dx = min(span * 0.035, 0.010)

    for pos in positions:
        if strand == "+":
            start, end = pos - dx, pos + dx
        elif strand == "-":
            start, end = pos + dx, pos - dx
        else:
            continue

        ax.add_patch(
            FancyArrowPatch(
                (start, y),
                (end, y),
                arrowstyle="-|>",
                mutation_scale=4.7,
                linewidth=0.45,
                color=GENE_LINE_COLOR,
                shrinkA=0,
                shrinkB=0,
                zorder=2,
            )
        )


def choose_shared_ymax(stage_data: dict, user_ymax: float | None = None) -> tuple[float, float]:
    all_neglog = pd.concat([d["NEGLOG10P"] for d in stage_data.values()]).astype(float)
    abs_max = float(all_neglog.max())

    if user_ymax is not None:
        ymin_needed = abs_max + max(0.25, abs_max * 0.03)
        return max(float(user_ymax), ymin_needed), abs_max

    # Always show all points: use true global maximum plus small headroom.
    pad = max(0.55, abs_max * 0.06)
    ymax = np.ceil((abs_max + pad) * 2.0) / 2.0

    return ymax, abs_max


def draw_clipped_points(ax, df: pd.DataFrame, ymax: float) -> None:
    """
    Draw clipped points as circles, not triangles.
    """
    clipped = df[df["NEGLOG10P"] > ymax].copy()

    if clipped.empty:
        return

    clipped_y = np.full(len(clipped), ymax - 0.16)
    has_r2 = clipped["r2"].notna()

    if (~has_r2).any():
        ax.scatter(
            clipped.loc[~has_r2, "POS"] / 1e6,
            clipped_y[~has_r2.to_numpy()],
            c=NO_LD_COLOR,
            s=POINT_SIZE,
            marker="o",
            linewidths=0,
            zorder=5,
        )

    if has_r2.any():
        ax.scatter(
            clipped.loc[has_r2, "POS"] / 1e6,
            clipped_y[has_r2.to_numpy()],
            c=clipped.loc[has_r2, "r2"].astype(float),
            cmap=LD_CMAP,
            norm=LD_NORM,
            s=POINT_SIZE,
            marker="o",
            linewidths=0,
            zorder=6,
        )


# ── Gene model handling ───────────────────────────────────────────────────────

def load_gene_models(
    gtf_path: str,
    chrom: str,
    start_bp: int,
    end_bp: int,
    max_genes: int = 30,
    show_uncharacterized: bool = False,
) -> list:
    open_fn = gzip.open if str(gtf_path).endswith(".gz") else open
    chr_set = _chr_aliases(chrom)
    genes = {}

    with open_fn(gtf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 9:
                continue

            seqname, _, feature, f_start, f_end, _, strand, _, attrs = parts

            if seqname not in chr_set:
                continue

            try:
                f_start_i = int(f_start)
                f_end_i = int(f_end)
            except ValueError:
                continue

            if f_end_i < start_bp or f_start_i > end_bp:
                continue

            at = parse_gtf_attr(attrs)
            gene_id = at.get("gene_id")

            if not gene_id:
                continue

            gene_name = at.get("gene_name", gene_id)

            if not should_plot_gene(
                gene_name,
                show_uncharacterized=show_uncharacterized,
            ):
                continue

            transcript_id = at.get("transcript_id")

            g = genes.setdefault(
                gene_id,
                {
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "strand": strand,
                    "start": f_start_i,
                    "end": f_end_i,
                    "transcripts": {},
                },
            )

            g["start"] = min(g["start"], f_start_i)
            g["end"] = max(g["end"], f_end_i)

            if transcript_id:
                t = g["transcripts"].setdefault(
                    transcript_id,
                    {
                        "start": f_start_i,
                        "end": f_end_i,
                        "exons": [],
                    },
                )

                t["start"] = min(t["start"], f_start_i)
                t["end"] = max(t["end"], f_end_i)

                if feature.lower() == "exon":
                    t["exons"].append((f_start_i, f_end_i))

    models = []

    for g in genes.values():
        if g["end"] < start_bp or g["start"] > end_bp:
            continue

        best_tx = None
        best_len = -1

        for tx in g["transcripts"].values():
            tx_len = tx["end"] - tx["start"] + 1

            if tx_len > best_len:
                best_len = tx_len
                best_tx = tx

        if best_tx:
            exons = sorted(best_tx["exons"]) if best_tx["exons"] else [(best_tx["start"], best_tx["end"])]
            m_start, m_end = best_tx["start"], best_tx["end"]
        else:
            exons = [(g["start"], g["end"])]
            m_start, m_end = g["start"], g["end"]

        models.append(
            {
                "gene_name": g["gene_name"],
                "strand": g["strand"],
                "start": m_start,
                "end": m_end,
                "exons": exons,
            }
        )

    models = sorted(models, key=lambda x: (x["start"], -(x["end"] - x["start"])))

    return models[:max_genes]


def assign_gene_lanes(models: list) -> list:
    lane_rightmost = []
    placed = []
    gap_bp = 50000

    for m in models:
        s, e = m["start"], m["end"]
        lane = None

        for i, right in enumerate(lane_rightmost):
            if s > right + gap_bp:
                lane = i
                lane_rightmost[i] = e
                break

        if lane is None:
            lane = len(lane_rightmost)
            lane_rightmost.append(e)

        mm = dict(m)
        mm["lane"] = lane
        placed.append(mm)

    return placed


# ── I/O ───────────────────────────────────────────────────────────────────────

def read_glm(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df.columns = [c.lstrip("#") for c in df.columns]

    if "TEST" in df.columns:
        df = df[df["TEST"] == "ADD"].copy()

    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df = df.dropna(subset=["POS", "P"])
    df = df[df["P"] > 0]
    df["NEGLOG10P"] = -np.log10(df["P"].astype(float))

    return df[["POS", "ID", "NEGLOG10P", "P"]].reset_index(drop=True)


def read_ld(path: str, target_loci: str) -> dict:
    is_plink1 = path.endswith(".ld")

    if is_plink1:
        df = pd.read_csv(
            path,
            sep=r"\s+",
            dtype=str,
            low_memory=False,
            skipinitialspace=True,
        )
        df.columns = [c.strip() for c in df.columns]
        id_a_col = "SNP_A"
        id_b_col = "SNP_B"
        r2_col = "R2"
    else:
        df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        df.columns = [c.lstrip("#") for c in df.columns]
        id_a_col = next((c for c in df.columns if c.upper() == "ID_A"), None)
        id_b_col = next((c for c in df.columns if c.upper() == "ID_B"), None)
        r2_col = next((c for c in df.columns if "R2" in c.upper()), None)

        if not all([id_a_col, id_b_col, r2_col]):
            sys.exit(
                f"ERROR: cannot locate ID_A/ID_B/R2 columns in {path}.\n"
                f"Columns present: {df.columns.tolist()}"
            )

    if r2_col not in df.columns or id_a_col not in df.columns or id_b_col not in df.columns:
        sys.exit(
            f"ERROR: expected columns {id_a_col}, {id_b_col}, {r2_col} in {path}.\n"
            f"Columns present: {df.columns.tolist()}"
        )

    target_set = _chr_variants(target_loci)
    ld = {}

    for _, row in df.iterrows():
        try:
            r2 = float(row.get(r2_col))
        except (ValueError, TypeError):
            continue

        if np.isnan(r2):
            continue

        a = str(row.get(id_a_col, ""))
        b = str(row.get(id_b_col, ""))

        if a in target_set:
            ld[b] = r2

        if b in target_set:
            ld[a] = r2

    for v in target_set:
        ld[v] = 1.0

    print(f"LD map: {len(ld)} variants loaded from {path}", file=sys.stderr)

    return ld


# ── Plot ──────────────────────────────────────────────────────────────────────

def plot_comparison(
    stage_data: dict,
    ld_maps: dict,
    target_loci: str,
    out_prefix: str,
    title: str,
    target_label: str = None,
    show_target_label: bool = True,
    genomic_locus: str = None,
    gene_gtf: str = None,
    max_track_genes: int = 30,
    show_uncharacterized_genes: bool = False,
    save_pdf: bool = False,
    y_max: float = None,
    fig_width: float = 6.25,
    fig_height: float = None,
) -> None:
    use_publication_style()

    target_label = target_label or target_loci

    stages = [s for s in STAGE_ORDER if s in stage_data]
    stages += sorted(s for s in stage_data if s not in STAGE_ORDER)
    n = len(stages)

    all_pos = pd.concat([d["POS"] for d in stage_data.values()])
    xmin, xmax = float(all_pos.min()), float(all_pos.max())

    # Smaller padding gives a compact regional x-range.
    pad = (xmax - xmin) * 0.006

    ymax, _ = choose_shared_ymax(stage_data, user_ymax=y_max)

    region = parse_genomic_locus(genomic_locus) if genomic_locus else None

    if region:
        region_chrom, region_start_bp, region_end_bp = region
    else:
        region_chrom = str(target_loci).split(":")[0] if ":" in str(target_loci) else "chr1"
        region_start_bp, region_end_bp = int(xmin), int(xmax)

    gene_models = []

    if gene_gtf:
        try:
            gene_models = load_gene_models(
                gtf_path=gene_gtf,
                chrom=region_chrom,
                start_bp=region_start_bp,
                end_bp=region_end_bp,
                max_genes=max_track_genes,
                show_uncharacterized=show_uncharacterized_genes,
            )
        except Exception as exc:
            print(f"WARNING: failed to parse gene GTF {gene_gtf}: {exc}", file=sys.stderr)
            gene_models = []

    nrows = n + (1 if gene_gtf else 0)

    # Association panels are intentionally taller; gene track remains compact.
    height_ratios = [1.12] * n + ([0.30] if gene_gtf else [])

    if fig_height is None:
        fig_height = 7.10 if gene_gtf else 6.40

    fig = plt.figure(figsize=(fig_width, fig_height), facecolor=FIGURE_FACE)

    # Balance panel compactness with clean separation.
    gs = GridSpec(
        nrows=nrows,
        ncols=1,
        figure=fig,
        height_ratios=height_ratios,
        hspace=0.12,
    )

    axes = []

    for idx in range(n):
        axes.append(fig.add_subplot(gs[idx, 0], sharex=axes[0] if axes else None))

    axes_all = list(axes)

    if gene_gtf:
        axes_all.append(fig.add_subplot(gs[-1, 0], sharex=axes[0]))

    genomewide_line = -np.log10(5e-8)

    try:
        target_pos_bp = int(target_loci.split(":")[1])
    except (IndexError, ValueError):
        target_pos_bp = None

    target_id_set = _chr_variants(target_loci)
    xlim = ((xmin - pad) / 1e6, (xmax + pad) / 1e6)

    major_locator = MaxNLocator(nbins=5, min_n_ticks=4, prune=None)

    for ax, stage in zip(axes, stages):
        df = stage_data[stage].copy()
        stage_ld = ld_maps.get(stage, {})
        df["r2"] = df["ID"].map(stage_ld)

        has_r2 = df["r2"].notna()
        no_r2 = ~has_r2

        if no_r2.any():
            ax.scatter(
                df.loc[no_r2, "POS"] / 1e6,
                df.loc[no_r2, "NEGLOG10P"],
                c=NO_LD_COLOR,
                s=POINT_SIZE,
                alpha=0.42,
                zorder=2,
                linewidths=0,
                marker="o",
                rasterized=True,
            )

        if has_r2.any():
            ax.scatter(
                df.loc[has_r2, "POS"] / 1e6,
                df.loc[has_r2, "NEGLOG10P"],
                c=df.loc[has_r2, "r2"].astype(float),
                cmap=LD_CMAP,
                norm=LD_NORM,
                s=POINT_SIZE,
                alpha=0.88,
                zorder=3,
                linewidths=0,
                marker="o",
                rasterized=True,
            )

        draw_clipped_points(ax, df, ymax)

        # Target variant: circular highlight.
        tgt_rows = df[df["ID"].isin(target_id_set)]

        if not tgt_rows.empty:
            ax.scatter(
                tgt_rows["POS"] / 1e6,
                tgt_rows["NEGLOG10P"],
                c=TARGET_COLOR,
                edgecolors="white",
                linewidths=0.40,
                s=POINT_SIZE,
                marker="D",
                zorder=6,
            )

        if target_pos_bp is not None:
            ax.axvline(
                target_pos_bp / 1e6,
                color=TARGET_GUIDE_COLOR,
                linewidth=0.56,
                linestyle=(0, (2.2, 2.8)),
                alpha=0.40,
                zorder=1,
            )

        ax.axhline(
            genomewide_line,
            color=THRESHOLD_COLOR,
            linewidth=0.60,
            linestyle=(0, (4, 3)),
            alpha=0.62,
        )

        ax.set_facecolor(PANEL_FACE)
        ax.grid(axis="y", color=GRID_COLOR, linewidth=0.52, linestyle="-", alpha=0.78)
        ax.grid(axis="x", visible=False)

        ax.set_ylabel(r"$-\log_{10}(\mathit{P})$")
        ax.set_ylim(0, ymax)
        ax.set_xlim(*xlim)

        ax.yaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=4))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.grid(axis="y", which="minor", color=GRID_COLOR, linewidth=0.40, linestyle="-", alpha=0.42)
        ax.xaxis.set_major_locator(major_locator)
        ax.xaxis.set_major_formatter(FuncFormatter(format_mb_ticks))

        ax.tick_params(axis="both", colors=TEXT_COLOR)
        ax.tick_params(axis="y", which="minor", length=2.0, width=0.45)
        ax.tick_params(axis="x", pad=2.8)

        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(True)
        ax.spines["left"].set_color(SPINE_COLOR)
        ax.spines["bottom"].set_color(SPINE_COLOR)
        ax.spines["top"].set_color(SPINE_COLOR)
        ax.spines["right"].set_color(SPINE_COLOR)
        ax.spines["left"].set_linewidth(0.70)
        ax.spines["bottom"].set_linewidth(0.70)
        ax.spines["top"].set_linewidth(0.70)
        ax.spines["right"].set_linewidth(0.70)

        ax.text(
            0.012,
            0.89,
            STAGE_LABELS.get(stage, stage),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9.1,
            fontweight="semibold",
            color=TEXT_COLOR,
            zorder=10,
        )

    for ax in axes[:-1]:
        ax.tick_params(axis="x", labelbottom=False)

    # When gene track exists, keep x tick labels only on the bottom-most gene axis.
    if gene_gtf and axes:
        axes[-1].tick_params(axis="x", labelbottom=False)

    if axes:
        fig.legend(
            handles=make_ld_legend_handles(),
            title=rf"LD $r^2$ with {target_label}",
            loc="upper left",
            bbox_to_anchor=(0.105, 0.962),
            frameon=False,
            borderaxespad=0.0,
            ncol=5,
            columnspacing=0.56,
            handlelength=0.60,
            handletextpad=0.22,
            labelspacing=0.16,
            title_fontsize=7.5,
            fontsize=7.4,
        )

        sig_handle = Line2D(
            [0],
            [0],
            color=THRESHOLD_COLOR,
            linewidth=0.58,
            linestyle=(0, (4, 3)),
            label=r"Genome-wide ($\mathit{P} < 5 \times 10^{-8}$)",
        )

        fig.legend(
            handles=[sig_handle],
            loc="upper right",
            bbox_to_anchor=(0.985, 0.962),
            frameon=False,
            borderaxespad=0.0,
            handlelength=1.35,
            handletextpad=0.28,
            fontsize=7.4,
        )

    if gene_gtf:
        gax = axes_all[-1]

        gax.set_facecolor("#FFFFFF")
        gax.set_xlim(*xlim)

        gax.spines["top"].set_visible(False)
        gax.spines["right"].set_visible(False)
        gax.spines["left"].set_visible(False)
        gax.spines["bottom"].set_color(SPINE_COLOR)
        gax.spines["bottom"].set_linewidth(0.68)

        gax.tick_params(axis="y", left=False, labelleft=False)
        gax.tick_params(axis="x", colors=TEXT_COLOR)
        gax.grid(False)

        gax.xaxis.set_major_locator(major_locator)
        gax.xaxis.set_major_formatter(FuncFormatter(format_mb_ticks))

        placed_genes = assign_gene_lanes(gene_models) if gene_models else []

        if placed_genes:
            lane_offsets = {}

            for gm in placed_genes:
                lane_y = gm["lane"] * 0.42
                line_y = lane_y + 0.18

                x0 = max(gm["start"], region_start_bp) / 1e6
                x1 = min(gm["end"], region_end_bp) / 1e6
                text_x = (x0 + x1) / 2.0

                gax.hlines(
                    y=line_y,
                    xmin=x0,
                    xmax=x1,
                    color=GENE_LINE_COLOR,
                    linewidth=0.60,
                    zorder=2,
                )

                draw_gene_direction_arrows(gax, x0, x1, line_y, gm["strand"])

                for es, ee in gm["exons"]:
                    ex0 = max(es, region_start_bp) / 1e6
                    ex1 = min(ee, region_end_bp) / 1e6

                    if ex1 <= ex0:
                        continue

                    gax.add_patch(
                        Rectangle(
                            (ex0, line_y - 0.055),
                            max(ex1 - ex0, 0.00085),
                            0.110,
                            facecolor=GENE_EXON_COLOR,
                            edgecolor="none",
                            alpha=0.95,
                            zorder=3,
                        )
                    )

                text_y = line_y - 0.10
                lane_key = gm["lane"]
                prev_end = lane_offsets.get(lane_key)

                if prev_end is not None and text_x - prev_end < 0.025:
                    text_x = prev_end + 0.025

                lane_offsets[lane_key] = text_x

                gax.text(
                    min(max(text_x, xlim[0] + 0.003), xlim[1] - 0.003),
                    text_y,
                    gm["gene_name"],
                    fontsize=7.4,
                    ha="center",
                    va="top",
                    color=TEXT_COLOR,
                    fontstyle="italic",
                    fontweight="medium",
                    clip_on=True,
                )

            gax.set_ylim(-0.16, max(gm["lane"] for gm in placed_genes) * 0.42 + 0.40)

        else:
            gax.set_ylim(0, 1)
            gax.text(
                (xlim[0] + xlim[1]) / 2,
                0.5,
                "No overlapping informative genes found in annotation",
                fontsize=8.0,
                ha="center",
                va="center",
                color=MUTED_TEXT_COLOR,
            )

        gax.set_xlabel(f"Genomic position on {region_chrom} (Mb)", labelpad=4.6)

    else:
        axes[-1].set_xlabel("Genomic position (Mb)", labelpad=5.0)

    axes_all[-1].set_xlim(*xlim)

    fig.suptitle(
        title,
        fontsize=11.0,
        fontweight="bold",
        x=0.5,
        y=0.983,
        ha="center",
        color=TEXT_COLOR,
    )

    fig.subplots_adjust(
        left=0.10,
        right=0.985,
        top=0.885,
        bottom=0.092,
    )

    out_path = out_prefix + ".png"
    plt.savefig(out_path, dpi=700, bbox_inches="tight")
    print(f"Saved: {out_path}")

    if save_pdf:
        pdf_path = out_prefix + ".pdf"
        plt.savefig(pdf_path, bbox_inches="tight")
        print(f"Saved: {pdf_path}")

    plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    stage_data = {}
    skipped = []

    for f in args.glm_files:
        stage = infer_stage(f)

        if stage is None:
            skipped.append(f)
            continue

        try:
            df = read_glm(f)
        except Exception as exc:
            print(f"WARNING: skipping {f}: {exc}", file=sys.stderr)
            continue

        if stage in stage_data:
            stage_data[stage] = pd.concat([stage_data[stage], df], ignore_index=True)
        else:
            stage_data[stage] = df

        print(f'Loaded {len(df)} variants for stage "{stage}" from {f}', file=sys.stderr)

    if skipped:
        print(f"WARNING: could not infer stage from: {skipped}", file=sys.stderr)

    if not stage_data:
        sys.exit("ERROR: no valid GLM files found.")

    ld_maps = {}

    for ld_path in args.ld_files:
        stage = infer_stage_from_ld(ld_path)

        if stage is None:
            print(
                f'WARNING: cannot infer stage from LD filename "{ld_path}"; '
                "expected ld_<stage>.vcor/.vcor1/.ld — skipping.",
                file=sys.stderr,
            )
            continue

        ld_maps[stage] = read_ld(ld_path, args.target_loci)

    missing_ld = [s for s in stage_data if s not in ld_maps]

    if missing_ld:
        print(
            f"WARNING: no LD file matched for stage(s): {missing_ld}. "
            "Those panels will be drawn entirely in grey.",
            file=sys.stderr,
        )

    plot_comparison(
        stage_data=stage_data,
        ld_maps=ld_maps,
        target_loci=args.target_loci,
        out_prefix=args.out,
        title=args.title,
        target_label=args.target_label,
        show_target_label=not args.hide_target_label,
        genomic_locus=args.genomic_locus,
        gene_gtf=args.gene_gtf,
        max_track_genes=args.max_track_genes,
        show_uncharacterized_genes=args.show_uncharacterized_genes,
        save_pdf=args.pdf,
        y_max=args.y_max,
        fig_width=args.fig_width,
        fig_height=args.fig_height,
    )


if __name__ == "__main__":
    main()
