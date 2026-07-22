#!/usr/bin/env python3
"""Genotype concordance + missingness of the down-sampled replicates vs the
cohort's own QC'd genotypes (refined_core), BEFORE and AFTER genotype QC.

Note on labels: the truth is NOT "30x" and the test is NOT "15x". Down-sampling is
per-platform, onto the cohort's baseline platform — the platforms it touches sit
anywhere from 18x to 35x, and each is scaled by its own measured fraction. So the
axes here are "cohort" (whatever depth each sample was sequenced at) and
"down-sampled" (all brought onto the baseline).

Produces ONE log and TWO figures:
  * confusion-matrix figure : 3 replicates x {pre_qc, post_qc} vs refined_core
  * missingness figure      : smiss / vmiss of the cohort vs ALL replicate
                              conditions (absolute + delta vs the cohort)

Confusion matrix (per the reference slides):
    rows = truth (refined_core, the cohort's own calls) in {RR, RA, AA}
    cols = test  (replicate, down-sampled)             in {RR, RA, AA, NC}
      RR = hom-ref (ALT dosage 0), RA = het (1), AA = hom-alt (2), NC = no-call
    (pairs where TRUTH is no-call are excluded and reported separately)
    The grid is every truth genotype, so a variant down-sampling deleted outright
    — gone from the replicate VCF, not merely no-called within it — lands in NC
    rather than escaping the tally.

Metrics (denominators follow the slides):
    called   = sum of all cells with test in {RR,RA,AA}
    NC       = sum of the NC column
    total    = called + NC
    Concordance rate    = (RR-RR + RA-RA + AA-AA) / called          (diagonal)
    Genotype miss rate  = NC / total
    False positive rate = (RR-RA + RR-AA + RA-AA) / called   (upper tri: test gained ALT)
    False negative rate = (RA-RR + AA-RR + AA-RA) / called   (lower tri: test lost ALT)
    (Concordance + FPR + FNR = 1 over the called cells)

Missingness (on the matrix grid — every truth variant, not only the ones each
replicate kept; a genotype absent from a dataset counts as missing):
    smiss[sample]  = missing genotypes of that sample / #variants
    vmiss[variant] = missing genotypes at that variant / #samples

Concordance and MissRate must be read together. Concordance is diagonal/called,
so it answers "when the replicate called, was it right" — genotypes it declined
to call are not in its denominator. Genotype QC pushes bad calls into NC, which
RAISES Concordance while raising MissRate. Neither number means much alone.

Inputs:
    --truth-raw  refined_core plink2 --export A .raw (ALT dosage 0/1/2/NA)
    --gt-dir     directory with per-condition long GT tables named
                 "<rep>.<cond>.gt.tsv" (cols: SAMPLE  ID  GT)
"""

import argparse
import glob
import os
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


GT_LABELS = ["RR", "RA", "AA", "NC"]
ORIG_LABEL = "cohort"
COND_ORDER = {"pre_qc": 0, "post_qc": 1}
COND_TITLE = {"pre_qc": "pre-GT-QC", "post_qc": "post-GT-QC"}


def fmt(n):
    """Integer with thousands separators."""
    return f"{int(n):,}"


def read_truth(path):
    """plink2 --export A .raw -> truth[iid][varid] = dosage in {0,1,2} or None."""
    truth = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        var_ids = [c.rsplit("_", 1)[0] for c in header[6:]]
        for row in fh:
            f = row.rstrip("\n").split("\t")
            iid = f[1]
            rec = {}
            for vid, val in zip(var_ids, f[6:]):
                rec[vid] = None if val == "NA" else int(round(float(val)))
            truth[iid] = rec
    return truth


def gt_to_dosage(gt):
    """VCF GT string -> ALT dosage {0,1,2} or None (no-call). Biallelic."""
    alleles = gt.replace("|", "/").split("/")
    if any(a == "." for a in alleles):
        return None
    return sum(1 for a in alleles if a != "0")


def read_test(path):
    """long GT table (SAMPLE ID GT) -> test[iid][varid] = dosage or None."""
    test = defaultdict(dict)
    with open(path) as fh:
        for row in fh:
            parts = row.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            test[parts[0]][parts[1]] = gt_to_dosage(parts[2])
    return test


def build_matrix(truth, test):
    """C[truth 0..2][test 0..3(NC)] over the truth grid; returns C, truth_nc, n.

    Every variant the truth carries is walked, not just the ones the replicate kept.
    Down-sampling can delete a variant outright — GenotypeGVCFs emits the forced
    site with ALT='.' when no sample supports the ALT any more, and normalization
    drops that row — and such a variant is simply absent from the replicate VCF.
    Absent means the replicate produced no genotype for it, i.e. a no-call for
    every sample, so it belongs in the NC column. Iterating over the intersection
    instead would drop it from the tally altogether and let a variant that
    vanished entirely cost nothing.

    (This does not move Concordance: NC is outside its denominator either way.
    It is MissRate that was understating the damage.)
    """
    C = np.zeros((3, 4), dtype=int)
    truth_nc = 0
    n = 0
    for iid in set(truth) & set(test):
        tr, xr = truth[iid], test[iid]
        for vid, t in tr.items():
            if t is None:
                truth_nc += 1
                continue
            x = xr.get(vid)          # absent from the replicate == no-call
            C[t][3 if x is None else x] += 1
            n += 1
    return C, truth_nc, n


def metrics(C):
    called = int(C[:, :3].sum())
    nc = int(C[:, 3].sum())
    total = called + nc
    diag = int(C[0, 0] + C[1, 1] + C[2, 2])
    upper = int(C[0, 1] + C[0, 2] + C[1, 2])   # test gained ALT -> FP
    lower = int(C[1, 0] + C[2, 0] + C[2, 1])   # test lost  ALT -> FN
    p = lambda a, b: (100.0 * a / b) if b else 0.0
    return dict(called=called, nc=nc, total=total, diag=diag, upper=upper, lower=lower,
                concordance=p(diag, called), miss_rate=p(nc, total),
                fpr=p(upper, called), fnr=p(lower, called))


def missingness(ds, samples, variants):
    """smiss (per sample) and vmiss (per variant) as % on the common grid.
    A genotype absent from the dataset counts as missing."""
    nv, ns = len(variants), len(samples)
    smiss, vmiss_cnt = [], defaultdict(int)
    for s in samples:
        rec = ds.get(s, {})
        m = 0
        for v in variants:
            if rec.get(v, None) is None:
                m += 1
                vmiss_cnt[v] += 1
        smiss.append(100.0 * m / nv if nv else 0.0)
    vmiss = [100.0 * vmiss_cnt[v] / ns if ns else 0.0 for v in variants]
    return np.array(smiss), np.array(vmiss)


def parse_label(fname):
    """'rep1.pre_qc.gt.tsv' -> (rep='rep1', cond='pre_qc')."""
    base = os.path.basename(fname)
    if base.endswith(".gt.tsv"):
        base = base[:-len(".gt.tsv")]
    rep, _, cond = base.partition(".")
    return rep, cond


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth-raw", required=True)
    ap.add_argument("--gt-dir", required=True)
    ap.add_argument("--out-log", required=True)
    ap.add_argument("--out-fig", required=True, help="confusion-matrix figure")
    ap.add_argument("--out-fig-missing", required=True, help="smiss/vmiss figure")
    args = ap.parse_args()

    truth = read_truth(args.truth_raw)
    tests = {}
    for tsv in sorted(glob.glob(os.path.join(args.gt_dir, "*.gt.tsv"))):
        tests[parse_label(tsv)] = read_test(tsv)

    results = {k: dict(zip(("C", "truth_nc", "n"), build_matrix(truth, t)))
               for k, t in tests.items()}
    for k in results:
        results[k]["m"] = metrics(results[k]["C"])

    reps = sorted({r for r, _ in results})
    conds = sorted({c for _, c in results}, key=lambda c: COND_ORDER.get(c, 9))
    order = [(r, c) for r in reps for c in conds if (r, c) in results]

    # ---- the grid missingness is measured on ----
    # Samples are intersected: only the down-sampled platforms are re-called, so a
    # sample no replicate carries has nothing to be missing from.
    samples = set(truth)
    for t in tests.values():
        samples &= set(t)
    samples = sorted(samples)
    # Variants are NOT intersected — the grid is the truth's, i.e. the matrix. A
    # variant a replicate lost outright is absent from that replicate, and absent
    # is exactly what "missing" means here (missingness() reads it via .get()).
    # Intersecting would delete such a variant from every dataset's grid, so the
    # one thing down-sampling did most violently would leave no trace.
    variants = sorted(set().union(*(set(rec) for rec in truth.values()))
                      if truth else set())

    miss = {ORIG_LABEL: missingness(truth, samples, variants)}
    for k, t in tests.items():
        miss[k] = missingness(t, samples, variants)

    # ---------------- combined log ----------------
    W = 96
    L = ["=" * W,
         " Genotype concordance & missingness vs refined_core (the cohort's own QC'd calls)",
         "=" * W,
         f" generated : {datetime.now().isoformat(timespec='seconds')}",
         f" truth     : {args.truth_raw}",
         " truth = refined_core (Array analogue); test = down-sampled replicate (WGS analogue)",
         f" shared grid: {fmt(len(samples))} samples x {fmt(len(variants))} variants"
         f" = {fmt(len(samples) * len(variants))} genotypes",
         "",
         " Metrics (over called cells unless noted; Concordance+FPR+FNR=1):",
         "   Concordance = diagonal / called      MissRate = NC / total",
         "   FPR = upper-tri (test gained ALT) / called   FNR = lower-tri (test lost ALT) / called",
         "",
         " CONCORDANCE SUMMARY",
         " " + "-" * (W - 2),
         f"   {'rep':<6}{'cond':<12}{'Concord%':>10}{'FPR%':>8}{'FNR%':>8}{'Miss%':>8}"
         f"{'called':>12}{'NC':>12}{'truthNC':>10}"]
    for rep, cond in order:
        r = results[(rep, cond)]
        m = r["m"]
        L.append(f"   {rep:<6}{COND_TITLE.get(cond, cond):<12}"
                 f"{m['concordance']:>10.3f}{m['fpr']:>8.3f}{m['fnr']:>8.3f}{m['miss_rate']:>8.3f}"
                 f"{fmt(m['called']):>12}{fmt(m['nc']):>12}{fmt(r['truth_nc']):>10}")

    L += ["",
          " MISSINGNESS SUMMARY  (mean over the common grid)",
          " " + "-" * (W - 2),
          f"   {'dataset':<20}{'smiss% mean':>13}{'smiss% max':>12}"
          f"{'vmiss% mean':>13}{'vmiss% max':>12}{'missing GT':>14}"]
    for key in [ORIG_LABEL] + order:
        name = key if key == ORIG_LABEL else f"{key[0]} {COND_TITLE.get(key[1], key[1])}"
        sm, vm = miss[key]
        n_missing = int(round(sm.mean() / 100.0 * len(variants) * len(samples)))
        L.append(f"   {name:<20}{sm.mean():>13.3f}{sm.max():>12.3f}"
                 f"{vm.mean():>13.3f}{vm.max():>12.3f}{fmt(n_missing):>14}")

    for rep, cond in order:
        r = results[(rep, cond)]
        C, m = r["C"], r["m"]
        L += [" " + "-" * (W - 2),
              f" {rep}  {COND_TITLE.get(cond, cond)}   Concordance={m['concordance']:.3f}%  "
              f"FPR={m['fpr']:.3f}%  FNR={m['fnr']:.3f}%  Miss={m['miss_rate']:.3f}%",
              "   confusion matrix (rows=truth cohort, cols=test down-sampled)",
              f"   {'':>9}{'RR':>13}{'RA':>13}{'AA':>13}{'NC':>13}"]
        for ti, tlab in enumerate(["RR", "RA", "AA"]):
            L.append(f"   {'truth ' + tlab:>9}" + "".join(f"{fmt(C[ti][xi]):>13}" for xi in range(4)))
        L.append(f"   truth no-call pairs excluded : {fmt(r['truth_nc'])}")
    L += ["=" * W, ""]
    Path(args.out_log).write_text("\n".join(L))
    sys.stdout.write("\n".join(L[:len(order) + 26]) + "\n")

    # ---------------- figure 1: confusion matrices ----------------
    nrow, ncol = len(reps), len(conds)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.0 * ncol, 3.9 * nrow), squeeze=False)
    for i, rep in enumerate(reps):
        for j, cond in enumerate(conds):
            ax = axes[i][j]
            r = results.get((rep, cond))
            if not r:
                ax.axis("off")
                continue
            C, m = r["C"], r["m"]
            rowsum = C.sum(axis=1, keepdims=True)
            frac = np.divide(C, rowsum, out=np.zeros_like(C, dtype=float), where=rowsum != 0)
            ax.imshow(frac, cmap="Oranges", vmin=0, vmax=1, aspect="auto")
            for ti in range(3):
                for xi in range(4):
                    ax.text(xi, ti, fmt(C[ti][xi]), ha="center", va="center", fontsize=9)
            ax.set_xticks(range(4)); ax.set_xticklabels(GT_LABELS)
            ax.set_yticks(range(3)); ax.set_yticklabels(["RR", "RA", "AA"])
            ax.set_xlabel("test (down-sampled)"); ax.set_ylabel("truth (cohort)")
            ax.set_title(f"{rep} · {COND_TITLE.get(cond, cond)}\n"
                         f"Conc {m['concordance']:.2f}%  Miss {m['miss_rate']:.2f}%\n"
                         f"FPR {m['fpr']:.3f}%  FNR {m['fnr']:.3f}%", fontsize=10)
            ax.axvline(2.5, color="grey", lw=1, ls="--")
    fig.suptitle("Down-sampled vs cohort (refined_core) genotype concordance",
                 fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(args.out_fig, dpi=150)
    plt.close(fig)

    # ---------------- figure 2: smiss / vmiss ----------------
    keys = [ORIG_LABEL] + order
    names = [ORIG_LABEL] + [f"{r}\n{COND_TITLE.get(c, c)}" for r, c in order]
    colors = ["0.6"] + ["#4C72B0" if c == "pre_qc" else "#C44E52" for _, c in order]

    fig2, ax2 = plt.subplots(2, 2, figsize=(14, 9))
    for row, (idx, mname) in enumerate([(0, "smiss (per sample)"), (1, "vmiss (per variant)")]):
        data = [miss[k][idx] for k in keys]
        bp = ax2[row][0].boxplot(data, labels=names, showfliers=False, patch_artist=True)
        for patch, col in zip(bp["boxes"], colors):
            patch.set_facecolor(col); patch.set_alpha(0.75)
        for i, d in enumerate(data):
            ax2[row][0].text(i + 1, d.mean(), f"{d.mean():.2f}", ha="center", va="bottom",
                             fontsize=8, fontweight="bold",
                             bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.0))
        ax2[row][0].set_ylabel(f"{mname}  missing %")
        ax2[row][0].set_title(f"{mname}: cohort vs all replicates")
        ax2[row][0].grid(axis="y", alpha=0.3)

        # delta vs original
        base = miss[ORIG_LABEL][idx]
        dd = [miss[k][idx] - base for k in order]
        bp2 = ax2[row][1].boxplot(dd, labels=names[1:], showfliers=False, patch_artist=True)
        for patch, col in zip(bp2["boxes"], colors[1:]):
            patch.set_facecolor(col); patch.set_alpha(0.75)
        for i, d in enumerate(dd):
            ax2[row][1].text(i + 1, d.mean(), f"{d.mean():+.2f}", ha="center", va="bottom",
                             fontsize=8, fontweight="bold",
                             bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.0))
        ax2[row][1].axhline(0, color="k", lw=1, ls="--")
        ax2[row][1].set_ylabel(f"Δ {mname}  (replicate − original) %")
        ax2[row][1].set_title(f"Δ {mname} vs cohort")
        ax2[row][1].grid(axis="y", alpha=0.3)

    fig2.suptitle(
        f"Genotype missingness: cohort (refined_core) vs down-sampled replicates\n"
        f"common grid {fmt(len(samples))} samples x {fmt(len(variants))} variants",
        fontsize=13)
    fig2.tight_layout(rect=[0, 0, 1, 0.95])
    fig2.savefig(args.out_fig_missing, dpi=150)
    plt.close(fig2)


if __name__ == "__main__":
    main()
