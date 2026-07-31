"""Shared publication style for the tuning.sample_qc figures.

One palette, one set of rcParams, a handful of helpers — so all figures read as one
system. Figures carry NO explanatory prose; each panel has a short title only. The
narrative lives in the figure captions in README.md (paper form).
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

# ---- palette (white ground, ink text, Okabe-Ito series) ----
SURFACE = "#ffffff"
INK     = "#1a1a1a"
MUTED   = "#5b5952"
GRID    = "#e7e5df"
AXIS    = "#9c988f"

C_CASE  = "#D55E00"   # vermillion
C_CTRL  = "#0072B2"   # blue
C_ACCENT = "#7B4FB5"  # purple (single-series highlights)
C_SIG   = "#B00020"   # threshold / warning red

# per-platform colours (HiSeqX = the control platform, kept near-black)
PCOL = {
    "HiSeqX 15x":        "#1a1a1a",
    "DNBSeq-T7 30x":     "#D55E00",
    "NovaSeq 30x":       "#009E73",
    "DNBSeq-G400RS 30x": "#0072B2",
    "DNBseq-G400RS 15x": "#CC79A7",
}
# soft tints for grouping bands
TINT_A = "#eef3f6"    # "accuracy" group
TINT_B = "#fbf1ea"    # "difference" group


def use_style():
    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor":   SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family":      "sans-serif",
        "font.sans-serif":  ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size":        10,
        "axes.titlesize":   11,
        "axes.titleweight": "regular",
        "axes.labelsize":   10,
        "axes.labelcolor":  MUTED,
        "xtick.labelsize":  9,
        "ytick.labelsize":  9,
        "xtick.color":      MUTED,
        "ytick.color":      MUTED,
        "axes.edgecolor":   AXIS,
        "axes.linewidth":   0.8,
        "axes.titlecolor":  INK,
        "legend.fontsize":  8.5,
        "legend.frameon":   False,
        "lines.solid_capstyle": "round",
        "figure.dpi":       120,
    })


def spine(ax, keep=("left", "bottom"), grid=None):
    """Despine to `keep`, colour ticks/spines, optional grid axis ('x'|'y'|'both')."""
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep)
    for s in keep:
        ax.spines[s].set_color(AXIS)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelcolor=MUTED, length=3)
    ax.set_axisbelow(True)
    if grid:
        ax.grid(axis=grid, color=GRID, linewidth=0.8)


def panel_title(ax, text, pad=6):
    ax.set_title(text, color=INK, fontsize=11, loc="left", pad=pad)


def panel_letter(ax, L, dx=-0.13, dy=1.03):
    """Bold panel letter outside the axes (paper convention)."""
    ax.text(dx, dy, L, transform=ax.transAxes, fontsize=13, fontweight="bold",
            color=INK, va="bottom", ha="left")


def save(fig, path, pad=0.30):
    fig.savefig(path, dpi=300, bbox_inches="tight", pad_inches=pad, facecolor=SURFACE)
    plt.close(fig)
