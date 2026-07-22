#!/usr/bin/env python3
"""What the down-sampling actually changed in the genotypes the association sees.

07_genotype_concordance scores the replicates on the 3,612-variant fixed_ready
matrix. That is the right set for asking "did re-calling at baseline depth
reproduce the cohort's genotypes", but it is not the set the association reads:
the SNP-based test reads random_model's variants in one region, the gene-based
test reads maf_lt_threshold's in another. This script scores the same comparison
on exactly those variants, so the number attaches to the actual result.

For each model x replicate it compares the association fileset against its base,
on the DOWN-SAMPLED samples only — every other sample is identical by
construction, and including them would dilute the disagreement with rows that
cannot disagree.

The confusion matrix and the metrics are imported from genotype_concordance.py
rather than restated, so the two reports cannot drift apart:
    rows = base (the cohort's own calls)      in {RR, RA, AA}
    cols = assoc (down-sampled swapped in)    in {RR, RA, AA, NC}
    Concordance = diagonal / called     MissRate = NC / total
    FPR = test gained ALT / called      FNR = test lost ALT / called
"""

import argparse
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from genotype_concordance import GT_LABELS, build_matrix, fmt, metrics, read_truth

W = 96


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", action="append", required=True, metavar="MODEL:REP:BASE_RAW:ASSOC_RAW",
                    help="one comparison; repeat. BASE_RAW/ASSOC_RAW are plink2 --export A .raw")
    ap.add_argument("--ds-samples", required=True, help="down-sampled sample IIDs, one per line")
    ap.add_argument("--out-log", required=True)
    ap.add_argument("--out-fig", required=True)
    args = ap.parse_args()

    ds = {l.strip() for l in open(args.ds_samples) if l.strip()}

    results = {}
    for spec in args.pair:
        model, rep, base_raw, assoc_raw = spec.split(":", 3)
        base = read_truth(base_raw)
        assoc = read_truth(assoc_raw)
        # Only the swapped samples: the rest are the same genotypes on both sides,
        # so they would pad the diagonal and flatter every metric.
        base = {s: v for s, v in base.items() if s in ds}
        assoc = {s: v for s, v in assoc.items() if s in ds}
        C, base_nc, n = build_matrix(base, assoc)
        results[(model, rep)] = dict(C=C, base_nc=base_nc, n=n, m=metrics(C),
                                     n_samples=len(set(base) & set(assoc)),
                                     n_vars=len(next(iter(base.values()))) if base else 0)

    models = sorted({m for m, _ in results})
    reps = sorted({r for _, r in results})
    order = [(m, r) for m in models for r in reps if (m, r) in results]

    L = ["=" * W,
         " What down-sampling changed in the association genotypes",
         "=" * W,
         f" generated : {datetime.now().isoformat(timespec='seconds')}",
         "",
         " Scored on the variants each association actually reads, and on the down-sampled",
         " samples only — every other sample carries its original genotypes untouched, so",
         " including them could only pad the diagonal.",
         "",
         "   rows = base  (the cohort's own calls)",
         "   cols = assoc (the replicate's re-called genotypes swapped in)",
         "   Concordance = diagonal / called      MissRate = NC / total",
         "   FPR = assoc gained ALT / called      FNR = assoc lost ALT / called",
         "",
         " SUMMARY",
         " " + "-" * (W - 2),
         f"   {'model':<12}{'rep':<6}{'samples':>9}{'variants':>10}{'Concord%':>10}"
         f"{'FPR%':>8}{'FNR%':>8}{'Miss%':>8}{'changed':>10}"]
    for model, rep in order:
        r = results[(model, rep)]
        m = r["m"]
        changed = m["called"] - m["diag"] + m["nc"]
        L.append(f"   {model:<12}{rep:<6}{r['n_samples']:>9,}{r['n_vars']:>10,}"
                 f"{m['concordance']:>10.3f}{m['fpr']:>8.3f}{m['fnr']:>8.3f}"
                 f"{m['miss_rate']:>8.3f}{fmt(changed):>10}")
    L += ["",
          "   'changed' counts every genotype that is not identical to the base: a different",
          "   call, or a call the replicate no longer makes. Those are the only cells that can",
          "   move an odds ratio.",
          ""]

    for model, rep in order:
        r = results[(model, rep)]
        C, m = r["C"], r["m"]
        L += [" " + "-" * (W - 2),
              f" {model}  {rep}   Concordance={m['concordance']:.3f}%  FPR={m['fpr']:.3f}%  "
              f"FNR={m['fnr']:.3f}%  Miss={m['miss_rate']:.3f}%",
              f"   {'':>9}{'RR':>13}{'RA':>13}{'AA':>13}{'NC':>13}"]
        for ti, tlab in enumerate(["RR", "RA", "AA"]):
            L.append(f"   {'base ' + tlab:>9}" + "".join(f"{fmt(C[ti][xi]):>13}" for xi in range(4)))
        L.append(f"   base no-call pairs excluded : {fmt(r['base_nc'])}")
    L += ["=" * W, ""]
    Path(args.out_log).write_text("\n".join(L))
    print("\n".join(L))

    # ── figure ──
    nrow, ncol = len(models), len(reps)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.8 * nrow), squeeze=False)
    for i, model in enumerate(models):
        for j, rep in enumerate(reps):
            ax = axes[i][j]
            r = results.get((model, rep))
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
            ax.set_xlabel("assoc (down-sampled)"); ax.set_ylabel("base (cohort)")
            ax.set_title(f"{model} · {rep}\n{r['n_samples']:,} samples x {r['n_vars']:,} variants\n"
                         f"Conc {m['concordance']:.2f}%  Miss {m['miss_rate']:.2f}%", fontsize=10)
            ax.axvline(2.5, color="grey", lw=1, ls="--")
    fig.suptitle("Association genotypes: down-sampled samples vs their original calls\n"
                 "on the variants each test actually reads",
                 fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(args.out_fig, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
