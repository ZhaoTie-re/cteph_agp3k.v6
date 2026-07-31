#!/usr/bin/env python3
"""Fig 4 — sample loss: current pipeline vs proposed within-platform design.
Clean panels only; narrative (incl. the rejected MAD variant) is in the caption.

OLD (as run): DP robust-Z within Target_DP < −3  OR  Het_F pooled mean ± 5·SD.
NEW (design): DP robust-Z within PLATFORM (MAD) < −3  OR  Het_F within-PLATFORM ± 5·SD.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from _style import (use_style, spine, panel_title, panel_letter, save,
                    INK, MUTED, C_CASE, C_CTRL, C_ACCENT)

MET = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/07_sample_qc/metrics/cteph_agp3k_v6_wgs_merged.sample_qc_metrics.tsv"
XLS = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/info/cteph_agp3k.v6.20260507.xlsx"
OUT = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tuning.sample_qc/figures/fig4_sample_loss.png"
C_GREY = "#b0ada5"

use_style()
m = pd.read_csv(MET, sep="\t", dtype={"IID": str})
info = pd.read_excel(XLS, dtype={"ID_JHRPv6": str})
info["ID_JHRPv6"] = info["ID_JHRPv6"].astype(str).str.strip()
j = m.merge(info[["ID_JHRPv6", "WGS_Platform", "Outcome"]], left_on="IID", right_on="ID_JHRPv6", how="left")
j["grp"] = np.where(j["Outcome"] == "PH", "case", "ctrl")
NCASE, NCTRL = (j.grp == "case").sum(), (j.grp == "ctrl").sum()
madz = lambda s: (s - s.median()) / (1.4826 * (s - s.median()).abs().median())
sdz = lambda s: (s - s.mean()) / s.std(ddof=0)
mu, sd = j["Het_F"].mean(), j["Het_F"].std()
old = (j.groupby("Target_DP")["DP"].transform(madz) < -3) | ((j["Het_F"] - mu).abs() > 5 * sd)
new = (j.groupby("WGS_Platform")["DP"].transform(madz) < -3) | (j.groupby("WGS_Platform")["Het_F"].transform(sdz).abs() > 5)
OLD, NEW = set(j.loc[old, "IID"]), set(j.loc[new, "IID"])
mad_reject = int((j.groupby("WGS_Platform")["Het_F"].transform(madz).abs() > 5).sum())
cc = lambda S: (int((j[j.IID.isin(S)].grp == "case").sum()), int((j[j.IID.isin(S)].grp == "ctrl").sum()))
old_c, old_k = cc(OLD); new_c, new_k = cc(NEW)
both, oo, no = OLD & NEW, OLD - NEW, NEW - OLD

fig, (a, b) = plt.subplots(1, 2, figsize=(11.6, 4.2), gridspec_kw=dict(width_ratios=[1.3, 1.0]))
fig.subplots_adjust(left=0.20, right=0.985, top=0.90, bottom=0.14, wspace=0.30)

# A — removed counts, stacked case/ctrl
rows = [("current pipeline\n(pooled Het, DP by target-depth)", old_c, old_k),
        ("proposed design\n(within-platform Het & DP)", new_c, new_k)]
for i, (lab, nc, nk) in enumerate(rows):
    a.barh(i, nk, color=C_CTRL, height=0.5); a.barh(i, nc, left=nk, color=C_CASE, height=0.5)
    a.text(nk + nc + 0.4, i, f"{nk+nc}  ({nc} case / {nk} ctrl)", va="center", fontsize=9, color=INK)
a.set_yticks([0, 1]); a.set_yticklabels([r[0] for r in rows]); a.set_ylim(-0.6, 1.6); a.invert_yaxis()
a.set_xlim(0, 34); a.set_xlabel("samples removed  (of 3,592)")
spine(a, keep=("bottom",), grid="x"); a.tick_params(axis="y", length=0, labelcolor=INK)
a.legend(handles=[Patch(color=C_CASE, label="case"), Patch(color=C_CTRL, label="control")], loc="lower right", labelcolor=INK)
a.text(0.985, 0.52, f"shared {len(both)}\nonly current {len(oo)} (case)\nonly proposed {len(no)} (ctrl)",
       transform=a.transAxes, fontsize=8.5, color=MUTED, va="center", ha="right", linespacing=1.5)
panel_title(a, f"Comparable total  ({len(OLD)} → {len(NEW)})"); panel_letter(a, "A", dx=-0.30)

# B — per-capita removal rate
xr = np.arange(2); w = 0.36
oldr = [100 * old_c / NCASE, 100 * old_k / NCTRL]; newr = [100 * new_c / NCASE, 100 * new_k / NCTRL]
b.bar(xr - w / 2, oldr, w, color=C_GREY); b.bar(xr + w / 2, newr, w, color=C_ACCENT)
for xx, v in zip(xr - w / 2, oldr): b.annotate(f"{v:.2f}%", (xx, v), xytext=(0, 3), textcoords="offset points", ha="center", fontsize=8, color=MUTED)
for xx, v in zip(xr + w / 2, newr): b.annotate(f"{v:.2f}%", (xx, v), xytext=(0, 3), textcoords="offset points", ha="center", fontsize=8, color=INK)
b.set_xticks(xr); b.set_xticklabels(["case", "control"]); b.set_ylabel("% of group removed")
spine(b, grid="y")
b.legend(handles=[Patch(color=C_GREY, label="current (pooled)"), Patch(color=C_ACCENT, label="proposed (within-platform)")], loc="upper right", labelcolor=INK)
panel_title(b, f"Less phenotype-differential  ({oldr[0]/oldr[1]:.1f}× → {newr[0]/newr[1]:.1f}×)"); panel_letter(b, "B", dx=-0.16)

save(fig, OUT)
print("wrote", OUT)
print(f"OLD={len(OLD)} ({old_c}c/{old_k}k)  NEW={len(NEW)} ({new_c}c/{new_k}k)  MAD-reject={mad_reject}")
