#!/usr/bin/env python3
"""Fig 2 — sample missingness (SMISS = 1 − call rate) is the confounded, benign axis.
Clean panels only; narrative is the figure caption in README.md."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from _style import (use_style, spine, panel_title, panel_letter, save,
                    INK, MUTED, SURFACE, C_CASE, C_CTRL, PCOL)

MET = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/07_sample_qc/metrics/cteph_agp3k_v6_wgs_merged.sample_qc_metrics.tsv"
XLS = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/info/cteph_agp3k.v6.20260507.xlsx"
OUT = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tuning.sample_qc/figures/fig2_smiss_confound.png"

use_style()
m = pd.read_csv(MET, sep="\t", dtype={"IID": str})
info = pd.read_excel(XLS, dtype={"ID_JHRPv6": str})
info["ID_JHRPv6"] = info["ID_JHRPv6"].astype(str).str.strip()
j = m.merge(info[["ID_JHRPv6", "WGS_Platform", "Outcome"]], left_on="IID", right_on="ID_JHRPv6", how="left")
j["grp"] = np.where(j["Outcome"] == "PH", "case", "ctrl")
j = j[j["SMISS"].notna()].copy()
order = j.groupby("WGS_Platform")["SMISS"].median().sort_values().index.tolist()
DEPTH = {"HiSeqX 15x": "15x", "DNBseq-G400RS 15x": "15x", "DNBSeq-T7 30x": "30x", "NovaSeq 30x": "30x", "DNBSeq-G400RS 30x": "30x"}

fig, (a, b) = plt.subplots(1, 2, figsize=(11.6, 4.4), gridspec_kw=dict(width_ratios=[1.1, 1.0]))
fig.subplots_adjust(left=0.16, right=0.985, top=0.90, bottom=0.13, wspace=0.28)

# A — SMISS box per platform
for i, p in enumerate(order):
    v = j.loc[j["WGS_Platform"] == p, "SMISS"].values; c = PCOL[p]
    a.boxplot([v], positions=[i], vert=False, widths=0.55, showfliers=False, patch_artist=True,
              boxprops=dict(facecolor=SURFACE, color=c, linewidth=1.8),
              whiskerprops=dict(color=c, linewidth=1.3), capprops=dict(color=c, linewidth=1.3),
              medianprops=dict(color=INK, linewidth=2.0))
a.set_yticks(range(len(order)))
a.set_yticklabels([f"{p}  [{DEPTH[p]}, {'ctrl' if p=='HiSeqX 15x' else 'case'}]" for p in order])
a.set_ylim(-0.6, len(order) - 0.4); a.invert_yaxis()
a.set_xlabel("SMISS  (sample missing rate = 1 − call rate)")
spine(a, keep=("bottom",), grid="x"); a.tick_params(axis="y", length=0, labelcolor=INK)
panel_title(a, "Missingness tracks depth"); panel_letter(a, "A", dx=-0.42)

# B — differential loss vs absolute threshold
thr = [0.03, 0.04, 0.05, 0.06, 0.08, 0.10]
nca, nct = (j.grp == "case").sum(), (j.grp == "ctrl").sum()
cp = [100 * ((j.SMISS > t) & (j.grp == "case")).sum() / nca for t in thr]
kp = [100 * ((j.SMISS > t) & (j.grp == "ctrl")).sum() / nct for t in thr]
x = np.arange(len(thr))
b.plot(x, kp, "-o", color=C_CTRL, lw=2, ms=6, label=f"control removed (n={nct:,})")
b.plot(x, cp, "-o", color=C_CASE, lw=2, ms=6, label=f"case removed (n={nca:,})")
for xi, c_, k_ in zip(x, cp, kp):
    b.annotate(f"{k_:.0f}%", (xi, k_), xytext=(0, 6), textcoords="offset points", ha="center", fontsize=7.5, color=C_CTRL)
    b.annotate(f"{c_:.0f}%", (xi, c_), xytext=(0, -12), textcoords="offset points", ha="center", fontsize=7.5, color=C_CASE)
b.set_xticks(x); b.set_xticklabels([f"{t:.2f}" for t in thr])
b.set_xlabel("absolute SMISS exclusion threshold"); b.set_ylabel("% of group removed")
spine(b, grid="both"); b.legend(loc="upper right", labelcolor=INK)
panel_title(b, "Any cut preferentially removes controls"); panel_letter(b, "B", dx=-0.16)

save(fig, OUT)
print("wrote", OUT)
print("SMISS median by platform:"); print(j.groupby("WGS_Platform")["SMISS"].median().round(3).to_string())
