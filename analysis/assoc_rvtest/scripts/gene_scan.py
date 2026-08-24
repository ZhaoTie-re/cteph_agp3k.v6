#!/usr/bin/env python3
# ==============================================================================
# gene_scan.py — scan QC and tier assignment for ONE gene-based association scan.
#
# The gene-level counterpart of _shared/scripts/call_peaks.py. One scan is one
# (impact stratum x test method) pair within one cohort.
#
# Why tiers come from FDR rather than a fixed P:
#   A gene-based scan tests a few hundred to a few thousand genes, and how many
#   depends on the impact stratum (164 genes under HIGH, ~8.6k under
#   MODERATE+HIGH, ~13k under LOW+MODERATE+HIGH in this study). A fixed P cut
#   would therefore mean a different multiple-testing burden in every stratum,
#   and the strata would stop being comparable. Benjamini-Hochberg adapts to the
#   number of genes actually tested, so the same FDR means the same expected
#   false-discovery proportion in every stratum.
#
# Outputs
#   gene_hits.tsv  one row per gene reaching a tier (the fan-out list)
#   gene_scan.tsv  every tested gene, with its tier — the full scan
#   scan_qc.tsv    one row: genes tested, lambda over gene P, tier counts
# ==============================================================================
import argparse
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

TIER_SIG, TIER_SUG, TIER_NONE = 'significant', 'suggestive', 'not_a_hit'

# THE PREFIX IS LOAD-BEARING. The pipeline derives a gene's tier from it
# (tierOf() in assoc_rvtest.nf), because a publishDir closure can only see
# process inputs and adding `tier` as an input would rehash every task. Changing
# these strings without changing tierOf() silently misfiles every result.
PREFIX = {TIER_SIG: 'sig', TIER_SUG: 'sug'}

# rvtest writes its gene-level statistics under these names; RANGE is
# 'chr:start-end' and is the only place the gene's position is given by some
# test methods.
P_COL = 'Pvalue'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--assoc', required=True, help='rvtest .assoc table (post NumVar filter)')
    p.add_argument('--cohort', required=True)
    p.add_argument('--stratum', required=True, help='impact stratum tag, e.g. moderate_high')
    p.add_argument('--method', required=True, help='test method tag, e.g. skato')
    p.add_argument('--fdr-significant', type=float, default=0.05)
    p.add_argument('--fdr-suggestive', type=float, default=0.10)
    p.add_argument('--out-dir', default='.')
    return p.parse_args()


def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values, NaN-safe and monotone."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    n = int(ok.sum())
    if not n:
        return out
    idx = np.flatnonzero(ok)
    order = idx[np.argsort(p[idx], kind='mergesort')]
    ranked = p[order] * n / np.arange(1, n + 1)
    # BH is the running minimum from the largest p downwards; without it the
    # adjusted values are not monotone and a gene can rank above one with a
    # smaller raw p.
    out[order] = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    return out


def lambda_gc(p):
    """Genomic-control inflation factor over the gene-level p-values."""
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    if p.size < 20:
        return np.nan
    return float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1))


def positions(df):
    """CHR/START/END, falling back to parsing RANGE when a method omits them."""
    have = {'CHR', 'START', 'END'} <= set(df.columns)
    if have and df[['CHR', 'START', 'END']].notna().all(axis=None):
        return df['CHR'].astype(str), df['START'], df['END']
    if 'RANGE' not in df.columns:
        raise SystemExit('assoc table has neither CHR/START/END nor RANGE')
    # 'chr4:76352000-76370000' or a comma-joined list of such spans; take the
    # outer bounds so a multi-span gene still gets one position.
    parsed = df['RANGE'].astype(str).str.extractall(
        r'(?P<c>[\w.]+):(?P<s>\d+)-(?P<e>\d+)')
    if not len(parsed):
        raise SystemExit('could not parse any span out of RANGE')
    g = parsed.groupby(level=0)
    return (g['c'].first().reindex(df.index),
            pd.to_numeric(g['s'].min()).reindex(df.index),
            pd.to_numeric(g['e'].max()).reindex(df.index))


def main():
    args = parse_args()
    df = pd.read_csv(args.assoc, sep='\t')
    if P_COL not in df.columns:
        raise SystemExit(f'{args.assoc}: no {P_COL} column (has: {list(df.columns)[:12]})')

    df[P_COL] = pd.to_numeric(df[P_COL], errors='coerce')
    n_read = len(df)
    # A gene whose test did not converge carries a blank or non-numeric p. It is
    # not a null result and must not enter the FDR denominator, so drop it here
    # and report the count rather than letting it dilute every other gene's
    # adjusted p.
    df = df[np.isfinite(df[P_COL]) & (df[P_COL] > 0) & (df[P_COL] <= 1)].copy()
    n_dropped = n_read - len(df)

    df['FDR_BH'] = bh_fdr(df[P_COL])
    # rvtest_post_process.py already writes an FDR column; recompute it here so
    # the tier decision has one owner, then check the two agree.
    if 'FDR' in df.columns:
        prior = pd.to_numeric(df['FDR'], errors='coerce')
        both = np.isfinite(prior) & np.isfinite(df['FDR_BH'])
        if both.any():
            gap = float(np.nanmax(np.abs(prior[both] - df['FDR_BH'][both])))
            if gap > 1e-6:
                print(f'[gene_scan] WARNING: incoming FDR differs from the recomputed BH '
                      f'value by up to {gap:.3g}; the tier uses the recomputed one',
                      file=sys.stderr)

    df['tier'] = np.select(
        [df['FDR_BH'] < args.fdr_significant, df['FDR_BH'] < args.fdr_suggestive],
        [TIER_SIG, TIER_SUG], default=TIER_NONE)

    chrom, start, end = positions(df)
    df['chrom'] = chrom.astype(str).str.replace('^chr', '', regex=True)
    df['start'] = pd.to_numeric(start, errors='coerce')
    df['end'] = pd.to_numeric(end, errors='coerce')
    df['cohort'] = args.cohort
    df['stratum'] = args.stratum
    df['method'] = args.method

    df = df.sort_values(P_COL, kind='mergesort').reset_index(drop=True)

    # gene_id: tier prefix + rank within the tier + a filesystem-safe gene name.
    # The rank makes the id stable and sortable; the name makes a directory
    # listing readable without a lookup.
    ids = []
    counters = {TIER_SIG: 0, TIER_SUG: 0}
    for _, r in df.iterrows():
        t = r['tier']
        if t == TIER_NONE:
            ids.append('')
            continue
        counters[t] += 1
        safe = re.sub(r'[^A-Za-z0-9_.-]', '_', str(r['Gene']))
        ids.append(f'{PREFIX[t]}{counters[t]:03d}_{safe}')
    df['gene_id'] = ids

    keep = ['cohort', 'stratum', 'method', 'gene_id', 'tier', 'Gene',
            'chrom', 'start', 'end', 'NumVar', 'NumPolyVar', 'N_INFORMATIVE',
            P_COL, 'FDR_BH']
    keep = [c for c in keep if c in df.columns]
    out = args.out_dir.rstrip('/')

    df[keep].to_csv(f'{out}/gene_scan.tsv', sep='\t', index=False)
    hits = df[df['tier'] != TIER_NONE]
    hits[keep].to_csv(f'{out}/gene_hits.tsv', sep='\t', index=False)

    lam = lambda_gc(df[P_COL])
    qc = pd.DataFrame([{
        'cohort': args.cohort, 'stratum': args.stratum, 'method': args.method,
        'n_genes_read': n_read, 'n_genes_dropped_no_p': n_dropped,
        'n_genes_tested': len(df),
        'lambda_gc': lam,
        'min_p': float(df[P_COL].min()) if len(df) else np.nan,
        'fdr_significant': args.fdr_significant, 'fdr_suggestive': args.fdr_suggestive,
        'n_significant': int((df['tier'] == TIER_SIG).sum()),
        'n_suggestive': int((df['tier'] == TIER_SUG).sum()),
        # The p at which BH crossed each tier, so the figure can draw the line
        # the decision actually used rather than a nominal one.
        'p_at_significant': _tier_edge(df, TIER_SIG),
        'p_at_suggestive': _tier_edge(df, [TIER_SIG, TIER_SUG]),
    }])
    qc.to_csv(f'{out}/scan_qc.tsv', sep='\t', index=False)

    print(f'[gene_scan] {args.cohort} {args.stratum}/{args.method}: '
          f'{len(df):,} genes tested ({n_dropped} dropped for a missing P), '
          f'lambda={lam:.3f}, {int((df["tier"] == TIER_SIG).sum())} significant, '
          f'{int((df["tier"] == TIER_SUG).sum())} suggestive')


def _tier_edge(df, tiers):
    """Largest raw p still inside the tier — the BH-implied threshold line."""
    tiers = [tiers] if isinstance(tiers, str) else tiers
    sel = df[df['tier'].isin(tiers)]
    return float(sel[P_COL].max()) if len(sel) else np.nan


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in gene_scan: {e}', file=sys.stderr)
        sys.exit(1)
