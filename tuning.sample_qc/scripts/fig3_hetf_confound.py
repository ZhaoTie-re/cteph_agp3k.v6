#!/usr/bin/env python3
"""Fig 3 — Het_F platform/phenotype confound, as distributions + η².
Clean panels only; narrative is the figure caption in README.md."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from _style import (use_style, spine, panel_title, panel_letter, save,
                    INK, MUTED, C_CASE, C_CTRL, C_ACCENT, PCOL)

MET = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/07_sample_qc/metrics/cteph_agp3k_v6_wgs_merged.sample_qc_metrics.tsv"
XLS = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/info/cteph_agp3k.v6.20260507.xlsx"
OUT = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tuning.sample_qc/figures/fig3_hetf_confound.png"

use_style()
m = pd.read_csv(MET, sep="\t", dtype={"IID": str})
info = pd.read_excel(XLS, dtype={"ID_JHRPv6": str})
info["ID_JHRPv6"] = info["ID_JHRPv6"].astype(str).str.strip()
j = m.merge(info[["ID_JHRPv6", "WGS_Platform", "Outcome"]], left_on="IID", right_on="ID_JHRPv6", how="left")
j["grp"] = np.where(j["Outcome"] == "PH", "case", "ctrl")
j = j[j["Het_F"].notna()].copy()
mu = j["Het_F"].mean()
order = j.groupby("WGS_Platform")["Het_F"].median().sort_values().index.tolist()

def eta2(group):
    g = j["Het_F"].mean(); sst = ((j["Het_F"] - g) ** 2).sum()
    ssb = sum(len(s) * (s["Het_F"].mean() - g) ** 2 for _, s in j.groupby(group))
    return ssb / sst
eta_plat, eta_grp = eta2("WGS_Platform"), eta2("grp")

xg = np.linspace(-0.06, 0.08, 600)
kde = lambda v: gaussian_kde(np.asarray(v))(xg)

fig = plt.figure(figsize=(11.6, 7.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.40, wspace=0.22,
                      left=0.07, right=0.985, top=0.94, bottom=0.08)

# A — per-platform densities
a = fig.add_subplot(gs[0, :])
a.axvline(mu, color=MUTED, lw=1.0, ls="--")
for p in order:
    d = kde(j.loc[j["WGS_Platform"] == p, "Het_F"]); n = int((j["WGS_Platform"] == p).sum())
    a.plot(xg, d, color=PCOL[p], lw=2.0, label=f"{p}  (n={n:,})"); a.fill_between(xg, d, color=PCOL[p], alpha=0.05)
a.set_xlim(-0.06, 0.08); a.set_xlabel("Het_F  (method-of-moments inbreeding coefficient)"); a.set_ylabel("density")
spine(a, grid="both"); panel_title(a, "Het_F density per platform"); panel_letter(a, "A", dx=-0.055)
a.legend(loc="upper right", labelcolor=INK)
a.text(mu, a.get_ylim()[1] * 0.97, " pooled mean 0.0086", color=MUTED, fontsize=8, va="top")

# B — case vs control densities
b = fig.add_subplot(gs[1, 0])
for g, c, lab in [("ctrl", C_CTRL, "control  (n=3,135)"),
                  ("case", C_CASE, "case  (n=457)")]:
    d = kde(j.loc[j["grp"] == g, "Het_F"]); b.plot(xg, d, color=c, lw=2.2, label=lab); b.fill_between(xg, d, color=c, alpha=0.08)
b.axvline(mu, color=MUTED, lw=1.0, ls="--")
b.set_xlim(-0.06, 0.08); b.set_xlabel("Het_F"); b.set_ylabel("density")
spine(b, grid="both"); panel_title(b, "Case vs control"); panel_letter(b, "B", dx=-0.13)
b.legend(loc="upper right", labelcolor=INK, handlelength=1.4, borderaxespad=0.4)

# C — variance explained (η²) with formula + meaning
c = fig.add_subplot(gs[1, 1])
for i, (lab, val, col) in enumerate([("by platform", eta_plat, C_ACCENT), ("by case / control", eta_grp, C_CASE)]):
    c.barh(i, 100 * val, color=col, height=0.5)
    c.text(100 * val + 0.05, i, f" {100*val:.2f}%", va="center", fontsize=11, color=INK, fontweight="bold")
c.set_yticks([0, 1]); c.set_yticklabels(["by platform", "by case / control"])
c.set_ylim(-0.6, 1.7); c.invert_yaxis(); c.set_xlim(0, 5); c.set_xlabel("% of Het_F variance explained")
spine(c, grid="x"); panel_title(c, "Variance explained (η²)"); panel_letter(c, "C", dx=-0.13)
c.text(0.985, 0.86,
       r"$\eta^2=\dfrac{SS_{between}}{SS_{total}}=\dfrac{\sum_g n_g(\bar x_g-\bar x)^2}{\sum_i (x_i-\bar x)^2}$",
       transform=c.transAxes, ha="right", va="top", fontsize=11, color=INK)
c.text(0.985, 0.58, "= share of Het_F spread BETWEEN groups\n(0 = groups identical). >5% = a batch effect.",
       transform=c.transAxes, ha="right", va="top", fontsize=8, color=MUTED, style="italic", linespacing=1.5)

save(fig, OUT)
print("wrote", OUT, f"| eta2 plat={100*eta_plat:.2f}% grp={100*eta_grp:.2f}%")
