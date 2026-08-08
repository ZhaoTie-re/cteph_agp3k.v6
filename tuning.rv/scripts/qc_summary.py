#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Per-minAC depth-confounding statistics on the BURDEN RATE
#           (rate = SMinAC / callable, per 1,000 callable sites):
#             - OLS on the rate in TWO declared models (one fit loop):
#                 noPC  rate ~ group + Z(MeanDP)        (primary, depth-adjusted)
#                 pc    rate ~ group + Z(MeanDP) + PCs  (+ancestry sensitivity)
#               giving beta_group (Case vs Control) and beta_depth.
#               SE = classical by default; --robust switches to HC3.
#             - Mann-Whitney U + Wasserstein of MeanDP across target depth.
#           Writes the one-row per-minAC stats TSV. (Poisson dropped: SMinAC is
#           a large near-normal sum and was massively over-dispersed.)
# Project : cteph_agp3k.v6 / tuning.rv  (rare-variant depth-confounding QC)
# Used by : tuning.rv.nf  process QC_SUMMARY
# ---------------------------------------------------------------------------
import argparse
import gzip
import subprocess
import sys

import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.formula.api as smf


def count_lines_fast(fs_path):
    """Count lines of a (gzipped) file via `gzip -dc | wc -l`, with a Python fallback."""
    try:
        with subprocess.Popen(['gzip', '-dc', fs_path], stdout=subprocess.PIPE) as p1:
            result = subprocess.check_output(['wc', '-l'], stdin=p1.stdout, text=True)
            if p1.stdout is not None:
                p1.stdout.close()
            p1.wait()
        return int(result.strip())
    except Exception as e:
        print(f"Warning: fast line count failed ({e}); slow fallback.", file=sys.stderr)
        count = 0
        opener = gzip.open if str(fs_path).endswith('.gz') else open
        with opener(fs_path, 'rt') as f:
            for i, _ in enumerate(f):
                count = i + 1
        return count


def load_pcs(pc_file, n_pcs):
    """Load ancestry PCs; return (DataFrame[SampleID, PC..], [pc_cols]). Empty if unavailable."""
    pc = pd.read_csv(pc_file, sep='\t', dtype=str)
    id_col = next((c for c in pc.columns if c.upper() in ('IID', '#IID', 'ID', 'SAMPLEID')), None)
    pc_cols = [c for c in pc.columns if c.lower().startswith('pc')][:n_pcs]
    if id_col is None or not pc_cols:
        raise ValueError(f"PC file lacks an ID column or pc columns: {list(pc.columns)}")
    for c in pc_cols:
        pc[c] = pd.to_numeric(pc[c], errors='coerce')
    pc = pc.rename(columns={id_col: 'SampleID'})
    pc['SampleID'] = pc['SampleID'].astype(str).str.strip()
    return pc[['SampleID'] + pc_cols], pc_cols


def _fit(formula, data, robust):
    """Fit OLS; return {'group':(b,l,u,p), 'depth':(b,l,u,p)}. NaN tuple on failure."""
    nan4 = (np.nan, np.nan, np.nan, np.nan)
    try:
        res = smf.ols(formula, data=data).fit(cov_type='HC3') if robust \
            else smf.ols(formula, data=data).fit()
        ci = res.conf_int()
        out = {}
        for term, key in (('Group_Bin', 'group'), ('MeanDP_Z', 'depth')):
            if term in res.params.index:
                out[key] = (res.params[term], ci.loc[term, 0], ci.loc[term, 1], res.pvalues[term])
            else:
                out[key] = nan4
        return out
    except Exception as e:
        print(f"Warning: OLS ({'HC3' if robust else 'classical'}) failed: {e}", file=sys.stderr)
        return {'group': nan4, 'depth': nan4}


def analyze_qc(args):
    try:
        variant_count = max(0, count_lines_fast(args.variant_metrics) - 1)

        df = pd.read_csv(args.sample_metrics, sep='\t', dtype={'SampleID': str})
        required = {'Group', 'MeanDP', 'SMinAC', 'TargetDP',
                    'SNumHomRef', 'SNumHet', 'SNumHomAlt'}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"sample metrics missing required columns: {sorted(missing)}")

        # ── Burden RATE (per 1,000 callable sites) ───────────────────────────
        df['callable'] = df['SNumHomRef'] + df['SNumHet'] + df['SNumHomAlt']
        df['rate_k'] = np.where(df['callable'] > 0, df['SMinAC'] / df['callable'] * 1000.0, np.nan)

        # Regression frame: Case/Control only, z-scored depth.
        reg = df[df['Group'].isin(['Case', 'Control'])].copy()
        reg['Group_Bin'] = reg['Group'].map({'Case': 1, 'Control': 0})
        sd = reg['MeanDP'].std()
        reg['MeanDP_Z'] = 0.0 if (sd == 0 or np.isnan(sd)) else (reg['MeanDP'] - reg['MeanDP'].mean()) / sd

        # (Spearman rho dropped: the depth-burden coupling is reported as the OLS
        #  depth coefficient beta_depth below — the single depth diagnostic.)

        # ── OLS models: declared once, fit via one loop ──────────────────────
        # noPC = primary (depth-adjusted); pc = +ancestry sensitivity (only when a
        # PC file is supplied). beta_group is emitted for every model; beta_depth
        # for the depth-adjusted ones.
        base = reg.dropna(subset=['rate_k', 'Group_Bin', 'MeanDP_Z'])
        nan_fit = {'group': (np.nan,) * 4, 'depth': (np.nan,) * 4}
        MODELS = [('noPC', True), ('pc', True)]        # (suffix, has_depth)
        formulas = {'noPC': "rate_k ~ Group_Bin + MeanDP_Z"}
        datasets = {'noPC': base}

        pc_source, n_pc_used = 'none', 0
        if args.pc_file:
            try:
                pc_df, pc_cols = load_pcs(args.pc_file, args.n_pcs)
                merged = base.merge(pc_df, on='SampleID', how='inner')
                if len(merged) > len(pc_cols) + 3:
                    formulas['pc'] = "rate_k ~ Group_Bin + MeanDP_Z + " + " + ".join(pc_cols)
                    datasets['pc'] = merged
                    pc_source, n_pc_used = args.pc_label, len(pc_cols)
            except Exception as e:
                print(f"Warning: PC model skipped ({e})", file=sys.stderr)

        fits = {sfx: (_fit(formulas[sfx], datasets[sfx], args.robust)
                      if (sfx in datasets and len(datasets[sfx]) > 2) else nan_fit)
                for sfx, _hd in MODELS}

        # ── Depth-distribution confounding across target depth ───────────────
        df['TargetDP_Str'] = df['TargetDP'].astype(str)
        groups = sorted(g for g in df['TargetDP_Str'].dropna().unique() if g.lower() != 'nan')
        mwu_p, wass = np.nan, np.nan
        if len(groups) >= 2:
            v1 = df.loc[df['TargetDP_Str'] == groups[0], 'MeanDP'].dropna().values
            v2 = df.loc[df['TargetDP_Str'] == groups[1], 'MeanDP'].dropna().values
            if len(v1) and len(v2):
                try:
                    _, mwu_p = stats.mannwhitneyu(v1, v2, alternative='two-sided')
                except Exception:
                    pass
                try:
                    wass = stats.wasserstein_distance(v1, v2)
                except Exception:
                    pass

        # ── Write one-row TSV (rate units: per 1,000 callable sites) ─────────
        row = {
            'MinAC_Threshold': args.min_ac,
            'Variant_Count': variant_count,
            'Sample_Count': len(df),
            'Rate_Mean': df['rate_k'].mean(),
        }
        for sfx, has_depth in MODELS:               # beta_group for all; beta_depth only for depth models
            g = fits[sfx]['group']
            (row[f'Beta_group_{sfx}'], row[f'Beta_group_{sfx}_L'],
             row[f'Beta_group_{sfx}_U'], row[f'P_group_{sfx}']) = g
            if has_depth:
                d = fits[sfx]['depth']
                (row[f'Beta_depth_{sfx}'], row[f'Beta_depth_{sfx}_L'],
                 row[f'Beta_depth_{sfx}_U'], row[f'P_depth_{sfx}']) = d
        row.update({
            'PC_source': pc_source, 'N_PC': n_pc_used,
            'SE_type': 'HC3' if args.robust else 'classical',
            'MWU_P_TargetDP': mwu_p, 'Wasserstein_Dist_TargetDP': wass,
        })
        cols = list(row.keys())
        with open(args.out_tsv, 'w') as f:
            f.write('\t'.join(cols) + '\n')
            f.write('\t'.join(f'{row[c]:.6g}' if isinstance(row[c], float) else str(row[c])
                              for c in cols) + '\n')
        print(f"Generated stats: {args.out_tsv}  (rate/1k; PC={pc_source}; SE={'HC3' if args.robust else 'classical'})")

    except Exception as e:
        print(f"Error in analyze_qc: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    p = argparse.ArgumentParser(description='Per-minAC depth-confounding stats on the burden rate.')
    p.add_argument('--sample-metrics', required=True)
    p.add_argument('--variant-metrics', required=True)
    p.add_argument('--min-ac', required=True, type=int)
    p.add_argument('--out-tsv', required=True)
    p.add_argument('--pc-file', default=None, help='Ancestry PC TSV (default bbj_mainland); enables the +PC model.')
    p.add_argument('--pc-label', default='bbj_mainland', help='Label recorded in PC_source.')
    p.add_argument('--n-pcs', type=int, default=10)
    p.add_argument('--robust', action='store_true', help='Use HC3 robust SE (default: classical).')
    args = p.parse_args()
    analyze_qc(args)


if __name__ == '__main__':
    main()
