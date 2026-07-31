#!/usr/bin/env python3
"""Sample QC filtering — within-platform design (see cteph_agp3k.v6/tuning.sample_qc/).

Two graded filters, both stratified WITHIN sequencing platform, plus a silent-pass fix:

  DP  (failed library)  : robust-Z (MAD) of Observed_Depth  < dp_threshold  (−3), one-sided low.
                          Precomputed per platform in the metrics table (DP_RobustZ_in_Platform).
  Het (excess-homozyg./ : Het_F within-platform  mean ± k·SD  (k=5, ddof=1), two-sided. The upper
       contamination)     tail is excess homozygosity; it is removed as "excess-homozygosity
                          outliers" (no ROH triage here).
  MISSING_METRIC        : a sample with no Het_F / DP / platform cannot be judged → QUARANTINED
                          (added to remove.id), never silently kept (fixes the old fillna(False)).

Sample-level call rate (SMISS) is NOT a filter here (it is the benign, depth/platform-confounded
axis); differential missingness is handled at the variant level. No sex check.

Output interface is unchanged: <prefix>.sample_qc.{detail.tsv, remove.id, keep.id, summary.json,
summary.txt, png}. remove.id / keep.id stay 2-col (#FID IID), no header.
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run within-platform sample QC from the metrics table.")
    parser.add_argument("--metrics-tsv", required=True, help="Input sample QC metrics TSV")
    parser.add_argument("--config-json", required=True, help="Sample QC config JSON")
    parser.add_argument("--out-prefix", required=True, help="Output file prefix")
    parser.add_argument(
        "--stratify-by",
        default="WGS_Platform",
        help="Metrics-table column used to stratify the DP and Het_F filters (default WGS_Platform).",
    )
    parser.add_argument("--sample-info-xlsx", help="Sample info Excel file for phenotype summary")
    parser.add_argument("--sample-id-col", help="Column name containing sample IDs in sample info")
    parser.add_argument("--phenotype-col", help="Column name containing phenotype labels in sample info")
    parser.add_argument("--case-value", default="PH", help="Case phenotype value")
    parser.add_argument("--ctrl-value", default="AGP3K", help="Control phenotype value")
    parser.add_argument("--case-label", help="Optional display label for cases; defaults to --case-value.")
    parser.add_argument("--ctrl-label", help="Optional display label for controls; defaults to --ctrl-value.")
    return parser.parse_args()


def normalize_text(value):
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    return text if text else None


def load_config(config_path: Path) -> dict:
    """New within-platform schema, with back-compat for the old condition1/condition2 keys.

    New form:
      { "stratify_by": "WGS_Platform",
        "dp_outlier": { "enabled": true, "robustz_threshold": -3.0 },
        "het_f":      { "enabled": true, "sd_multiplier": 5.0 },
        "quarantine_missing_metrics": true }
    """
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = json.load(handle)
    if not isinstance(cfg, dict):
        cfg = {}

    dp = cfg.get("dp_outlier", {}) if isinstance(cfg.get("dp_outlier"), dict) else {}
    het = cfg.get("het_f", {}) if isinstance(cfg.get("het_f"), dict) else {}
    # Back-compat: old condition1.dp_robustz.remove_if_less_than / condition2.sd_multiplier
    old_c1 = cfg.get("condition1", {}) if isinstance(cfg.get("condition1"), dict) else {}
    old_c2 = cfg.get("condition2", {}) if isinstance(cfg.get("condition2"), dict) else {}

    dp_threshold = dp.get(
        "robustz_threshold",
        old_c1.get("dp_robustz", {}).get("remove_if_less_than", cfg.get("dp_robustz_threshold", -3.0)),
    )
    het_k = het.get(
        "sd_multiplier",
        old_c2.get("sd_multiplier", cfg.get("het_f_sd_multiplier", 5.0)),
    )
    return {
        "stratify_by": cfg.get("stratify_by", "WGS_Platform"),
        "dp_enabled": bool(dp.get("enabled", True)),
        "dp_threshold": float(dp_threshold),
        "het_enabled": bool(het.get("enabled", old_c2.get("enabled", True))),
        "het_sd_multiplier": float(het_k),
        "quarantine_missing_metrics": bool(cfg.get("quarantine_missing_metrics", True)),
    }


def read_metrics(metrics_path: Path, stratify_by: str) -> pd.DataFrame:
    df = pd.read_csv(metrics_path, sep="\t", dtype=str)

    # DP robust-Z column: prefer the per-platform column; accept the legacy name for safety.
    rz_col = None
    for cand in ("DP_RobustZ_in_Platform", "DP_RobustZ_in_TargetDP"):
        if cand in df.columns:
            rz_col = cand
            break
    if rz_col is None:
        raise ValueError("Metrics table lacks a DP robust-Z column (DP_RobustZ_in_Platform).")
    if rz_col != "DP_RobustZ_in_Platform":
        df["DP_RobustZ_in_Platform"] = df[rz_col]

    required = ["#FID", "IID", "Het_F", "DP", "SMISS", "DP_RobustZ_in_Platform", stratify_by]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in metrics table: {', '.join(missing)}")

    for c in ("Het_F", "DP", "SMISS", "DP_RobustZ_in_Platform"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df[stratify_by] = df[stratify_by].map(normalize_text)
    if "Target_DP" not in df.columns:
        df["Target_DP"] = None
    return df


def attach_phenotype(df, xlsx_path, sample_id_col, phenotype_col, case_value, ctrl_value):
    """Attach PHENO_LABEL / PHENO_GROUP (case/ctrl/other/unknown) from the sample sheet."""
    result = df.copy()
    result["PHENO_LABEL"] = None
    result["PHENO_GROUP"] = "unknown"
    if not xlsx_path or not xlsx_path.exists() or not sample_id_col or not phenotype_col:
        return result

    info = pd.read_excel(xlsx_path, engine="openpyxl")
    if sample_id_col not in info.columns or phenotype_col not in info.columns:
        return result

    sub = info[[sample_id_col, phenotype_col]].copy()
    sub.columns = ["IID", "PHENO_RAW"]
    sub["IID"] = sub["IID"].map(normalize_text)
    sub["PHENO_RAW"] = sub["PHENO_RAW"].map(normalize_text)
    sub = sub[sub["IID"].notna()].copy()

    dup = sub["IID"].duplicated(keep=False)
    if dup.any():
        conflict = sub.loc[dup].drop_duplicates(subset=["IID", "PHENO_RAW"], keep=False)["IID"].unique().tolist()
        if conflict:
            print(f"WARNING: conflicting phenotype for {len(conflict)} IID(s); ignored in summary.", file=sys.stderr)
        sub = sub[~sub["IID"].isin(conflict)]

    result["PHENO_LABEL"] = result["IID"].map(sub.set_index("IID")["PHENO_RAW"])
    case_norm, ctrl_norm = normalize_text(case_value), normalize_text(ctrl_value)

    def _group(label):
        if label is None or (isinstance(label, float) and pd.isna(label)):
            return "unknown"
        if label == case_norm:
            return "case"
        if label == ctrl_norm:
            return "ctrl"
        return "other"

    result["PHENO_GROUP"] = result["PHENO_LABEL"].map(_group)
    return result


def apply_filters(df: pd.DataFrame, config: dict) -> tuple:
    """Within-platform DP (one-sided low) + Het_F (two-sided) + quarantine of missing metrics."""
    strat = config["stratify_by"]
    k = config["het_sd_multiplier"]
    out = df.copy()

    # (0) Quarantine: a sample with no Het_F / DP / platform cannot be judged. These are removed
    #     with reason MISSING_METRIC (never silently kept), and excluded from the per-platform stats.
    missing = out["Het_F"].isna() | out["DP"].isna() | out[strat].isna()
    if not config["quarantine_missing_metrics"]:
        missing = pd.Series(False, index=out.index)
    out["QC_FAIL_MISSING"] = missing
    judged = ~missing

    # (1) DP outlier: per-platform robust-Z already in the table; one-sided low.
    if config["dp_enabled"]:
        fail_dp = judged & (out["DP_RobustZ_in_Platform"] < config["dp_threshold"])
    else:
        fail_dp = pd.Series(False, index=out.index)
    out["QC_FAIL_DP"] = fail_dp.fillna(False)

    # (2) Het_F outlier: per-platform mean ± k·SD (ddof=1), two-sided.
    out["HET_PLAT_MEAN"] = np.nan
    out["HET_PLAT_LOW"] = np.nan
    out["HET_PLAT_HIGH"] = np.nan
    plat_bounds = {}
    if config["het_enabled"]:
        for plat, idx in out.loc[judged].groupby(strat).groups.items():
            vals = out.loc[idx, "Het_F"]
            mean = vals.mean(skipna=True)
            sd = vals.std(skipna=True, ddof=1)
            if pd.isna(mean) or pd.isna(sd) or sd == 0:
                lo = hi = np.nan
            else:
                lo, hi = mean - k * sd, mean + k * sd
            plat_bounds[plat] = (mean, sd, lo, hi)
            out.loc[idx, "HET_PLAT_MEAN"] = mean
            out.loc[idx, "HET_PLAT_LOW"] = lo
            out.loc[idx, "HET_PLAT_HIGH"] = hi
        fail_het = judged & (
            (out["Het_F"] < out["HET_PLAT_LOW"]) | (out["Het_F"] > out["HET_PLAT_HIGH"])
        )
    else:
        fail_het = pd.Series(False, index=out.index)
    out["QC_FAIL_HET"] = fail_het.fillna(False)

    out["QC_REMOVE"] = out["QC_FAIL_MISSING"] | out["QC_FAIL_DP"] | out["QC_FAIL_HET"]

    def _reason(row):
        r = []
        if row["QC_FAIL_MISSING"]:
            r.append("MISSING_METRIC")
        if row["QC_FAIL_DP"]:
            r.append("DP_OUTLIER")
        if row["QC_FAIL_HET"]:
            r.append("HET_OUTLIER")
        return ";".join(r)

    out["QC_REASON"] = out.apply(_reason, axis=1)

    stats = {
        "design": "within_platform",
        "stratify_by": strat,
        "n_total": int(len(out)),
        "n_remove": int(out["QC_REMOVE"].sum()),
        "n_keep": int((~out["QC_REMOVE"]).sum()),
        "n_fail_dp": int(out["QC_FAIL_DP"].sum()),
        "n_fail_het": int(out["QC_FAIL_HET"].sum()),
        "n_quarantine_missing": int(out["QC_FAIL_MISSING"].sum()),
        "dp_enabled": config["dp_enabled"],
        "dp_threshold": config["dp_threshold"],
        "het_enabled": config["het_enabled"],
        "het_sd_multiplier": k,
        "per_platform_het": {
            str(p): {
                "n": int((out[strat] == p).sum()),
                "mean": None if pd.isna(m) else float(m),
                "sd": None if pd.isna(s) else float(s),
                "low": None if pd.isna(lo) else float(lo),
                "high": None if pd.isna(hi) else float(hi),
            }
            for p, (m, s, lo, hi) in plat_bounds.items()
        },
    }
    return out, stats


def write_outputs(df: pd.DataFrame, stats: dict, out_prefix: Path) -> None:
    detail_path = Path(f"{out_prefix}.sample_qc.detail.tsv")
    remove_path = Path(f"{out_prefix}.sample_qc.remove.id")
    keep_path = Path(f"{out_prefix}.sample_qc.keep.id")
    summary_json_path = Path(f"{out_prefix}.sample_qc.summary.json")
    summary_txt_path = Path(f"{out_prefix}.sample_qc.summary.txt")

    df.to_csv(detail_path, sep="\t", index=False)

    # 2-col, no header — consumed by `plink2 --remove`. Interface must stay identical.
    df.loc[df["QC_REMOVE"], ["#FID", "IID"]].to_csv(remove_path, sep="\t", header=False, index=False)
    df.loc[~df["QC_REMOVE"], ["#FID", "IID"]].to_csv(keep_path, sep="\t", header=False, index=False)

    with summary_json_path.open("w", encoding="utf-8") as handle:
        json.dump(stats, handle, ensure_ascii=False, indent=2)

    with summary_txt_path.open("w", encoding="utf-8") as handle:
        for key in ("n_total", "n_remove", "n_keep", "n_fail_dp", "n_fail_het", "n_quarantine_missing"):
            handle.write(f"{key}\t{stats[key]}\n")
        handle.write(f"stratify_by\t{stats['stratify_by']}\n")
        handle.write(f"dp_threshold\t{stats['dp_threshold']}\n")
        handle.write(f"het_sd_multiplier\t{stats['het_sd_multiplier']}\n")


# ----------------------------------------------------------------------------------
# Visualization — publication house style from cteph_agp3k.v6/tuning.sample_qc/.
# ----------------------------------------------------------------------------------
SURFACE, INK, MUTED, GRID, AXIS = "#ffffff", "#1a1a1a", "#5b5952", "#e7e5df", "#9c988f"
C_CASE, C_CTRL, C_KEEP = "#C8102E", "#005A9C", "#cfcdc6"   # case / control / retained-bar grey
C_OTHER, C_SIG = "#8a877f", "#8a877f"
BAND, C_CLOUD = "#f3f2ee", "#dcdad3"   # neutral mean±SD band; faint kept-point cloud


def _pheno_color(group: str) -> str:
    return {"case": C_CASE, "ctrl": C_CTRL}.get(group, C_OTHER)


def _spine(ax, keep=("left", "bottom"), grid=None):
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep)
    for s in keep:
        ax.spines[s].set_color(AXIS)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=MUTED, length=3)
    ax.set_axisbelow(True)
    if grid:
        ax.grid(axis=grid, color=GRID, linewidth=0.8)


def generate_visualization(df, stats, case_value, ctrl_value, out_prefix, case_label=None, ctrl_label=None):
    """Publication-grade 4-panel figure:
      A  Het_F within platform (per-platform mean ± k·SD window); removed points coloured case/ctrl.
      B  DP robust-Z within platform (threshold line); removed points coloured case/ctrl.
      C  Samples removed, by reason × phenotype (grouped bars).
      D  Retained vs removed, per phenotype (stacked bars with counts and %).
    """
    from matplotlib.lines import Line2D

    strat = stats["stratify_by"]
    k = stats["het_sd_multiplier"]
    thr = stats["dp_threshold"]
    case_label = str(case_label if case_label is not None else case_value)
    ctrl_label = str(ctrl_label if ctrl_label is not None else ctrl_value)

    full = df.copy()
    if "PHENO_GROUP" not in full.columns:
        full["PHENO_GROUP"] = "unknown"
    d = full[full[strat].notna()].copy()  # scatter panels need a platform axis
    if d.empty:
        fig, ax = plt.subplots(figsize=(8, 3), facecolor=SURFACE)
        ax.text(0.5, 0.5, "no platform-labelled samples to plot", ha="center", va="center", color=MUTED)
        ax.axis("off")
        plt.savefig(Path(f"{out_prefix}.sample_qc.png"), dpi=300, bbox_inches="tight", facecolor=SURFACE)
        plt.close()
        return

    order = d.groupby(strat)["Het_F"].median().sort_values().index.tolist()
    ypos = {p: i for i, p in enumerate(order)}
    rng = np.random.default_rng(0)

    # ---- case/ctrl counts, from the FULL table (includes any platform-missing quarantined) ----
    def cc(mask):
        s = full[mask]
        return int((s["PHENO_GROUP"] == "case").sum()), int((s["PHENO_GROUP"] == "ctrl").sum())
    n_case = int((full["PHENO_GROUP"] == "case").sum())
    n_ctrl = int((full["PHENO_GROUP"] == "ctrl").sum())
    rem_case, rem_ctrl = cc(full["QC_REMOVE"])
    keep_case, keep_ctrl = n_case - rem_case, n_ctrl - rem_ctrl
    reasons = [("Het outlier", full["QC_FAIL_HET"]),
               ("DP outlier", full["QC_FAIL_DP"]),
               ("Missing metric", full["QC_FAIL_MISSING"])]

    fig = plt.figure(figsize=(13.6, 4.0 + 1.05 * len(order)), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05 * len(order) + 0.4, 2.1],
                          hspace=0.40, wspace=0.22, left=0.165, right=0.975, top=0.905, bottom=0.065)
    a, b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    c, e = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

    from scipy.stats import gaussian_kde

    C_RAIN = "#8a857b"          # retained points — warm grey, cohesive with the cloud
    CLOUD_F, CLOUD_E, C_MED = "#d7d3c9", "#a09c90", "#6d6a61"   # density fill / edge / median

    def _raincloud(ax, sub, xcol, i, rem_col):
        """Raincloud per row: half-violin density (cloud, with a median line) above the baseline,
        jittered retained points (rain) below, removed points coloured by phenotype."""
        vals = pd.to_numeric(sub[xcol], errors="coerce").dropna().values
        if len(vals) > 2 and np.nanstd(vals) > 1e-9:
            kde = gaussian_kde(vals)
            xs = np.linspace(vals.min(), vals.max(), 256)
            raw = kde(xs)
            dmax = float(raw.max())
            dd = raw / dmax * 0.42
            ax.fill_between(xs, i + 0.08, i + 0.08 + dd, color=CLOUD_F, edgecolor=CLOUD_E,
                            linewidth=1.1, zorder=2)
            med = float(np.median(vals))
            ax.plot([med, med], [i + 0.08, i + 0.08 + float(kde(med)[0]) / dmax * 0.42],
                    color=C_MED, lw=1.3, zorder=3)
        kept, rem = sub[~sub["QC_REMOVE"]], sub[sub[rem_col]]
        ax.scatter(kept[xcol], i - 0.14 - rng.uniform(0.0, 0.20, len(kept)),
                   s=5.5, color=C_RAIN, alpha=0.45, linewidths=0, zorder=4)
        ax.scatter(rem[xcol], i - 0.24 + rng.uniform(-0.03, 0.03, len(rem)), s=60,
                   color=[_pheno_color(g) for g in rem["PHENO_GROUP"]], edgecolors="white",
                   linewidths=1.2, zorder=6)

    ylo, yhi = -0.55, len(order) - 1 + 0.66

    # ---- A: Het_F raincloud per platform; warm band = per-platform mean ± k·SD window ----
    for p in order:
        i = ypos[p]
        sub = d[d[strat] == p]
        lo = sub["HET_PLAT_LOW"].dropna().unique()
        hi = sub["HET_PLAT_HIGH"].dropna().unique()
        if len(lo) and len(hi) and np.isfinite(lo[0]) and np.isfinite(hi[0]):
            a.add_patch(plt.Rectangle((lo[0], i - 0.38), hi[0] - lo[0], 0.92, facecolor="#f6f2ea",
                                       edgecolor="#e7ddc9", linewidth=0.9, zorder=0))
        _raincloud(a, sub, "Het_F", i, "QC_FAIL_HET")
    a.set_yticks(range(len(order)))
    a.set_yticklabels([f"{p}  (n={int((d[strat]==p).sum()):,})" for p in order], fontsize=9.5)
    a.set_ylim(ylo, yhi)
    a.set_xlabel("heterozygosity  $F$", color=INK, fontsize=11)
    _spine(a, keep=("bottom",), grid="x")
    a.tick_params(axis="x", labelsize=10)
    a.tick_params(axis="y", length=0, labelcolor=INK)
    a.set_title("A", color=INK, fontsize=13, loc="left", fontweight="bold", pad=8)
    a.set_title(f"Heterozygosity by platform  ·  band = mean ± {k:g} SD", color=INK, fontsize=11, loc="center", pad=8)

    # ---- B: depth robust-Z raincloud per platform, with a clearly-marked removal threshold ----
    for p in order:
        _raincloud(b, d[d[strat] == p], "DP_RobustZ_in_Platform", ypos[p], "QC_FAIL_DP")
    zvals = pd.to_numeric(d["DP_RobustZ_in_Platform"], errors="coerce")
    xmin, xmax = min(float(zvals.min()), thr) - 0.7, float(zvals.max()) + 0.4
    b.set_xlim(xmin, xmax)
    b.axvspan(xmin, thr, color="#f6e3e6", alpha=0.85, zorder=0)                 # removal zone
    b.axvline(thr, color=C_SIG, lw=1.8, ls=(0, (6, 2)), zorder=6)
    b.set_yticks(range(len(order)))
    b.set_yticklabels([])
    b.set_ylim(ylo, yhi)
    b.set_xlabel("sequencing-depth robust-$Z$", color=INK, fontsize=11)
    _spine(b, keep=("bottom",), grid="x")
    b.tick_params(axis="x", labelsize=10)
    b.tick_params(axis="y", length=0)
    ticks = sorted(set([t for t in range(-2, int(xmax) + 1, 2)] + [int(thr)]))
    b.set_xticks(ticks)
    for lbl, tv in zip(b.get_xticklabels(), ticks):
        if abs(tv - thr) < 1e-9:
            lbl.set_color(C_SIG)
            lbl.set_fontweight("bold")
    b.set_title("B", color=INK, fontsize=13, loc="left", fontweight="bold", pad=8)
    b.set_title(f"Depth outlier by platform  ·  remove if $Z < {thr:g}$", color=INK, fontsize=11, loc="center", pad=8)
    b.annotate(f"$Z<{thr:g}$: remove", xy=(thr, yhi), xytext=(-4, -3), textcoords="offset points",
               color=C_SIG, fontsize=8.5, ha="right", va="top", fontweight="bold")

    # ---- C: removed by reason × phenotype ----
    labels = [r[0] for r in reasons]
    casec = [cc(m)[0] for _, m in reasons]
    ctrlc = [cc(m)[1] for _, m in reasons]
    y, hh = np.arange(len(labels)), 0.36
    c.barh(y - hh / 2, casec, hh, color=C_CASE, label=f"case ({case_label})")
    c.barh(y + hh / 2, ctrlc, hh, color=C_CTRL, label=f"control ({ctrl_label})")
    for yy, v in list(zip(y - hh / 2, casec)) + list(zip(y + hh / 2, ctrlc)):
        c.text(v + max(casec + ctrlc + [1]) * 0.02, yy, f"{v}", va="center", fontsize=8.5,
               color=INK if v else MUTED)
    c.set_yticks(y)
    c.set_yticklabels(labels, fontsize=9.5)
    c.set_ylim(-0.6, len(labels) - 0.4)
    c.invert_yaxis()
    c.set_xlim(0, max(casec + ctrlc + [1]) * 1.28)
    c.set_xlabel("samples removed", color=MUTED, fontsize=10)
    _spine(c, keep=("bottom", "left"), grid="x")
    c.tick_params(axis="y", length=0, labelcolor=INK)
    c.set_title("C", color=INK, fontsize=12, loc="left", fontweight="bold", pad=6)
    c.set_title("Removed, by reason", color=INK, fontsize=10.5, loc="center", pad=6)

    # ---- D: retained vs removed, per phenotype (each bar normalised to its own group total, so
    #         the 10×-smaller case group is still readable) ----
    rows = [("control", n_ctrl, keep_ctrl, rem_ctrl, C_CTRL),
            ("case", n_case, keep_case, rem_case, C_CASE)]
    for i, (lab, tot, kp, rm, col) in enumerate(rows):
        if tot == 0:
            continue
        kpf, rmf = 100.0 * kp / tot, 100.0 * rm / tot
        e.barh(i, kpf, color=C_KEEP, height=0.55, zorder=2)
        e.barh(i, rmf, left=kpf, color=col, height=0.55, zorder=3)
        e.text(kpf * 0.5, i, f"retained {kp:,}  ({kpf:.1f}%)", va="center", ha="center", fontsize=8.5, color=INK)
        e.text(101, i, f"removed {rm}", va="center", ha="left", fontsize=9, color=col, fontweight="bold")
    e.set_yticks([0, 1])
    e.set_yticklabels([f"control (n={n_ctrl:,})", f"case (n={n_case:,})"], fontsize=9.5)
    e.set_ylim(-0.6, 1.6)
    e.set_xlim(0, 120)
    e.set_xticks([0, 25, 50, 75, 100])
    e.set_xlabel("% of group", color=MUTED, fontsize=10)
    _spine(e, keep=("bottom",), grid="x")
    e.tick_params(axis="y", length=0, labelcolor=INK)
    e.set_title("D", color=INK, fontsize=12, loc="left", fontweight="bold", pad=6)
    e.set_title("Retained vs removed", color=INK, fontsize=10.5, loc="center", pad=6)

    fig.legend(handles=[
        Line2D([], [], marker="o", ls="None", color="#7f7c73", markersize=7, label="retained"),
        Line2D([], [], marker="o", ls="None", color=C_CASE, markersize=8, markeredgecolor="white", label="removed (case)"),
        Line2D([], [], marker="o", ls="None", color=C_CTRL, markersize=8, markeredgecolor="white", label="removed (control)"),
    ], loc="upper right", frameon=False, fontsize=10, ncol=3, bbox_to_anchor=(0.975, 0.997), labelcolor=INK,
        handletextpad=0.4, columnspacing=1.5)
    fig.suptitle(
        f"Sample-level QC  ·  within-platform stratification  ·  {stats['n_remove']} of {stats['n_total']:,} removed",
        color=INK, fontsize=12.5, fontweight="bold", x=0.02, ha="left", y=0.978)
    plt.savefig(Path(f"{out_prefix}.sample_qc.png"), dpi=300, bbox_inches="tight", facecolor=SURFACE)
    plt.close()


def main() -> int:
    args = parse_args()
    try:
        config = load_config(Path(args.config_json))
        strat = args.stratify_by or config["stratify_by"]
        config["stratify_by"] = strat
        metrics = read_metrics(Path(args.metrics_tsv), strat)
        if args.sample_info_xlsx and args.sample_id_col and args.phenotype_col:
            metrics = attach_phenotype(
                metrics, Path(args.sample_info_xlsx), args.sample_id_col, args.phenotype_col,
                args.case_value, args.ctrl_value,
            )
        filtered, stats = apply_filters(metrics, config)
        write_outputs(filtered, stats, Path(args.out_prefix))
        case_label = args.case_label or args.case_value
        ctrl_label = args.ctrl_label or args.ctrl_value
        generate_visualization(filtered, stats, args.case_value, args.ctrl_value,
                               Path(args.out_prefix), case_label, ctrl_label)
        print(
            f"Sample QC (within-platform) finished: total={stats['n_total']}, remove={stats['n_remove']} "
            f"(DP={stats['n_fail_dp']}, Het={stats['n_fail_het']}, missing={stats['n_quarantine_missing']}), "
            f"keep={stats['n_keep']}",
            file=sys.stderr,
        )
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
