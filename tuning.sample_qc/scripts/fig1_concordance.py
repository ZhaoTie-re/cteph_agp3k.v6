#!/usr/bin/env python3
"""Fig 1 — WGS-vs-array concordance by depth stratum (T7-30x vs HiSeqX-15x).
Clean panels only; narrative is the figure caption in README.md."""
import os, sys, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from _style import (use_style, spine, panel_letter, save, INK, MUTED, C_CTRL, C_CASE,
                    TINT_A, TINT_B, C_SIG)

PKL = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/tuning.concordance.v5/results/05.merge_concordance_vmiss_summary/merged._DP8_GQ20_LAF0.2_HAF0.8_.summary.pkl"
OUT = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tuning.sample_qc/figures/fig1_concordance.png"

use_style()
obj = pickle.load(open(PKL, "rb"))
row = {[x for x in k if x in ("ALL", "15X", "30X")][0]: v[3].iloc[0] for k, v in obj.items()}

# (letter, title, column, ylim, fmt, is_difference-panel)
M = [("A", "Genotype concordance", "GENOTYPE_CONCORDANCE", (0.995, 1.0),   "{:.4f}", False),
     ("B", "False-positive rate",  "FALSE_POSITIVE_RATE",  (0, 2.2e-4),    "{:.1e}", False),
     ("C", "False-negative rate",  "FALSE_NEGATIVE_RATE",  (0, 2.2e-4),    "{:.1e}", False),
     ("D", "Genotype miss rate",   "GENOTYPE_MISS_RATE",   (0, 0.05),      "{:.3f}", True)]

fig, axes = plt.subplots(1, 4, figsize=(11.6, 3.5))
fig.subplots_adjust(left=0.05, right=0.99, top=0.86, bottom=0.16, wspace=0.42)
for ax, (L, title, col, ylim, fmt, diff) in zip(axes, M):
    ax.set_facecolor(TINT_B if diff else TINT_A)
    v15, v30 = row["15X"][col], row["30X"][col]
    ax.bar(0, v15, 0.60, color=C_CTRL)
    ax.bar(1, v30, 0.60, color=C_CASE)
    for x, val in [(0, v15), (1, v30)]:
        ax.annotate(fmt.format(val), (x, val), xytext=(0, 3), textcoords="offset points",
                    ha="center", fontsize=9, color=INK, fontweight="bold")
    ax.set_ylim(*ylim); ax.set_xlim(-0.6, 1.6)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["HiSeqX\n15x", "T7\n30x"], fontsize=8.5)
    spine(ax, grid="y")
    ax.set_title(title, color=INK, fontsize=10.5, loc="left", pad=6)
    panel_letter(ax, L, dx=-0.30, dy=1.05)
axes[3].text(0.96, 0.96, "1.8× more\nat 15x", transform=axes[3].transAxes, ha="right",
             va="top", fontsize=10, color=C_SIG, fontweight="bold", linespacing=1.25)
save(fig, OUT)
print("wrote", OUT)
for s in ("15X", "30X"):
    r = row[s]
    print(f"{s}: conc={r['GENOTYPE_CONCORDANCE']:.5f} FP={r['FALSE_POSITIVE_RATE']:.2e} "
          f"FN={r['FALSE_NEGATIVE_RATE']:.2e} miss={r['GENOTYPE_MISS_RATE']:.4f}")
