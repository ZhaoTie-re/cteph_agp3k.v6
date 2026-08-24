#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Quantify the relatedness the GRM encodes, for one cohort.
#
#           THE POINT OF THE MIXED MODEL IS THAT RELATED SAMPLES ARE KEPT. The
#           fixed-effects analysis removes them; this one retains them and lets
#           the GRM account for the correlation. That trade is the single biggest
#           design difference between the two components, and until it is
#           measured it is only an assertion — this makes it a number.
#
#           Computed on EXACTLY the markers the GRM is built from (the LD-pruned
#           set from PRUNE_MARKERS), so it describes the matrix SAIGE actually
#           fitted rather than a differently-ascertained proxy.
#
#           THE METRIC IS THE GRM ITSELF — the variance-standardized relationship
#           matrix plink2 computes with --make-rel, which is the same quantity
#           SAIGE builds from the same markers. So these are the values the model
#           actually used, not a proxy for them.
#
#           Two parts are reported separately because they mean different things:
#             OFF-DIAGONAL  relatedness between two samples. Degree bands are the
#                           usual expectations: 0.5 (1st), 0.25 (2nd), 0.125 (3rd),
#                           ~1.0 (duplicate/MZ) — twice the KING kinship scale,
#                           because a GRM relationship is 2 x kinship.
#             DIAGONAL      a sample's relatedness to itself, i.e. 1 + inbreeding.
#                           SAIGE OVERRIDES THIS with --isDiagofKinSetAsOne=True,
#                           so what it was before that is provenance the fit
#                           otherwise discards.
#
#           CAVEAT, stated rather than avoided: a GRM value conflates relatedness
#           with shared ancestry, so in ancestry-filtered cohorts the off-diagonal
#           bulk shifts with the filter and is not a pure kinship estimate. The
#           degree bands are therefore indicative, not diagnostic.
# Component: assoc_saige
# Used by : assoc_saige.nf  process RELATEDNESS
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Expected GRM relationship by degree, strongest first; a pair takes the first
# class it clears. These are 2 x the KING kinship thresholds because a GRM
# relationship is twice a kinship coefficient.
DEGREES = [('mz_or_dup', 0.708, 'duplicate / MZ twin'),
           ('first',     0.354, '1st degree (parent-offspring, full sib)'),
           ('second',    0.177, '2nd degree (half sib, grandparent, avuncular)'),
           ('third',     0.0884, '3rd degree (first cousin)')]
UNRELATED_BELOW = 0.0884
# Histogram of the FULL pairwise distribution, so the figure never has to carry
# the matrix itself.
# Uniform 0.005 resolution across the whole range. The tail must be binned as
# finely as the bulk: coarse bins above the third-degree boundary turn the
# informative part of the distribution into a few flat steps.
HIST_EDGES = np.arange(-0.60, 1.205, 0.01)
# The diagonal sits near 1 and needs its own, finer, scale.
DIAG_EDGES = np.arange(0.5, 1.505, 0.005)


def parse_args():
    p = argparse.ArgumentParser(description='Relatedness among one cohort, on the GRM markers.')
    p.add_argument('--rel-bin', required=True, help='plink2 --make-rel triangle bin4 output')
    p.add_argument('--rel-id', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--fam', required=True, help='the analysed (random-model) .fam')
    p.add_argument('--compare-fam', help='the fixed-effects .fam, to count what relatedness cost it')
    p.add_argument('--min-relationship', type=float, default=UNRELATED_BELOW,
                   help='report pairs at or above this GRM relationship (default: 3rd degree).')
    p.add_argument('--out-pairs', required=True)
    p.add_argument('--out-summary', required=True)
    p.add_argument('--out-hist', required=True)
    p.add_argument('--out-diag', required=True)
    return p.parse_args()


def read_triangle(path, n):
    """plink2's lower triangle, float32 -> (off_diagonal, diagonal).

    `--make-rel triangle` INCLUDES the diagonal (row i holds (i,0)..(i,i)), while
    `--make-king triangle` excludes it. Both sizes are accepted and told apart by
    the element count, so the file is never silently misread as the other.
    """
    v = np.fromfile(path, dtype=np.float32)
    with_diag = n * (n + 1) // 2
    without = n * (n - 1) // 2
    if v.size == with_diag:
        # Row i holds i+1 entries, (i,0)..(i,i), so row i starts at i(i+1)/2 and
        # its diagonal element (i,i) is at i(i+1)/2 + i = i(i+3)/2.
        i = np.arange(n, dtype=np.int64)
        idx = i * (i + 3) // 2
        diag = v[idx]
        # A variance-standardized diagonal is 1 + inbreeding, so it sits near 1.
        # If the index were wrong it would pick off-diagonal values, which centre
        # on 0 — cheap to check and it turns a silent mis-read into an error.
        med = float(np.median(diag))
        if not (0.5 < med < 1.5):
            raise SystemExit(f'{path}: extracted diagonal has median {med:.4f}, which is not a '
                             f'relationship-with-self — the triangle indexing is wrong')
        off = np.delete(v, idx)
        return off, diag
    if v.size == without:
        return v, None
    raise SystemExit(f'{path}: {v.size:,} values; {n:,} samples need '
                     f'{with_diag:,} (with diagonal) or {without:,} (without)')


def pair_indices(n):
    """Row/column index of every entry of the lower triangle, in file order."""
    r = np.repeat(np.arange(1, n), np.arange(1, n))
    c = np.concatenate([np.arange(i) for i in range(1, n)])
    return r, c


def main():
    args = parse_args()
    ids = pd.read_csv(args.rel_id, sep='\t')
    ids.columns = [c.lstrip('#') for c in ids.columns]
    idcol = 'IID' if 'IID' in ids.columns else ids.columns[-1]
    sid = ids[idcol].astype(str).to_numpy()
    n = len(sid)

    k, diag = read_triangle(args.rel_bin, n)

    # ── the off-diagonal distribution, as a histogram ──────────────────────
    counts, edges = np.histogram(k, bins=HIST_EDGES)
    pd.DataFrame({'cohort': args.cohort, 'lo': edges[:-1], 'hi': edges[1:],
                  'n_pairs': counts}).to_csv(args.out_hist, sep='\t', index=False)

    # ── the diagonal: self-relatedness, before SAIGE forces it to 1 ────────
    if diag is not None:
        dc, de = np.histogram(diag, bins=DIAG_EDGES)
        pd.DataFrame({'cohort': args.cohort, 'lo': de[:-1], 'hi': de[1:],
                      'n_samples': dc}).to_csv(args.out_diag, sep='\t', index=False)
    else:
        pd.DataFrame(columns=['cohort', 'lo', 'hi', 'n_samples']).to_csv(
            args.out_diag, sep='\t', index=False)

    # ── the pairs that matter ─────────────────────────────────────────────
    sel = np.flatnonzero(k >= args.min_relationship)
    r, c = pair_indices(n)
    kin = k[sel]
    # np.select, not a loop of assignments: DEGREES is ordered strongest-first
    # and every threshold is implied by the ones above it, so successive
    # `deg[kin >= lo] = name` assignments each overwrite the previous and
    # everything ends up in the weakest class.
    deg = np.select([kin >= lo for _n, lo, _l in DEGREES],
                    [n for n, _lo, _l in DEGREES], default='unrelated')
    pairs = pd.DataFrame({'cohort': args.cohort, 'id1': sid[r[sel]], 'id2': sid[c[sel]],
                          'grm': kin, 'degree': deg}).sort_values('grm', ascending=False)
    pairs.to_csv(args.out_pairs, sep='\t', index=False)

    # ── who is involved, and what it would have cost to drop them ─────────
    involved = set(pairs['id1']) | set(pairs['id2'])
    fam = pd.read_csv(args.fam, sep=r'\s+', header=None,
                      names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'], dtype={1: str})
    fam['IID'] = fam['IID'].astype(str)
    n_case = int((fam['PHENO'] == 2).sum())
    inv_case = int(fam[fam['IID'].isin(involved)]['PHENO'].eq(2).sum())

    row = {'cohort': args.cohort, 'n_samples': n, 'n_cases': n_case,
           'n_controls': int((fam['PHENO'] == 1).sum()),
           'n_pairs_total': int(n * (n - 1) // 2),
           'grm_offdiag_median': float(np.median(k)), 'grm_offdiag_max': float(k.max()),
           'grm_diag_median': float(np.median(diag)) if diag is not None else np.nan,
           'grm_diag_min': float(diag.min()) if diag is not None else np.nan,
           'grm_diag_max': float(diag.max()) if diag is not None else np.nan,
           'grm_diag_note': 'computed; SAIGE refits with --isDiagofKinSetAsOne=True',
           'n_related_pairs': int(len(pairs)),
           'n_samples_with_relative': len(involved),
           'pct_samples_with_relative': 100.0 * len(involved) / n if n else np.nan,
           'n_cases_with_relative': inv_case}
    for name, _lo, _label in DEGREES:
        row[f'n_pairs_{name}'] = int((pairs['degree'] == name).sum())

    # The direct statement of what the GRM buys: the samples the fixed-effects
    # component had to drop, which this one keeps.
    if args.compare_fam and Path(args.compare_fam).exists():
        other = pd.read_csv(args.compare_fam, sep=r'\s+', header=None,
                            names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'], dtype={1: str})
        kept = set(other['IID'].astype(str))
        dropped = set(fam['IID']) - kept
        row['n_absent_from_fixed_model'] = len(dropped)
        row['n_absent_and_related'] = len(dropped & involved)
        row['n_cases_absent_from_fixed_model'] = int(
            fam[fam['IID'].isin(dropped)]['PHENO'].eq(2).sum())

    pd.DataFrame([row]).to_csv(args.out_summary, sep='\t', index=False)
    dtxt = (f'; diagonal median {np.median(diag):.4f} [{diag.min():.4f}, {diag.max():.4f}]'
            if diag is not None else '')
    print(f'[relatedness] {args.cohort}: {n:,} samples, {len(pairs):,} pairs at GRM '
          f'>= {args.min_relationship:g} involving {len(involved):,} samples '
          f'({row["pct_samples_with_relative"]:.1f}%); max off-diagonal {k.max():.4f}{dtxt}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in relatedness: {e}', file=sys.stderr)
        sys.exit(1)
