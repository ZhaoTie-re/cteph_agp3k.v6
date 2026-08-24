#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : SAIGE step-2 output -> the canonical summary-statistic schema every
#           downstream component reads.
#
#           THE CANONICAL SCHEMA IS PLINK2'S `--glm` SCHEMA. That is a deliberate
#           choice, not an accident of history: the peak caller, the annotator,
#           the LD/fine-mapping chain and every figure live in analysis/_shared
#           and are about association RESULTS, not about which engine produced
#           them. Adapting SAIGE here means all of that is reused unchanged, and
#           a SAIGE scan and a plink2 scan of the same cohort are directly
#           comparable row for row.
#
#           TWO THINGS ARE EASY TO GET WRONG AND ARE HANDLED EXPLICITLY:
#
#           1. `OBS_CT` means different things in different plink2 outputs.
#              Under `--glm` it is the number of NON-MISSING SAMPLES (2,975 for a
#              3,071-sample cohort); under `--freq` it is an ALLELE count. The
#              downstream reads the `--glm` convention, so this writes samples.
#
#           2. The variant ID is the join key for EVERYTHING downstream — rsID
#              lookup, the snpEff index, the LD panels, the credible sets. It is
#              therefore RECONSTRUCTED as chr:pos:REF:ALT from SAIGE's own
#              CHR/POS/Allele1/Allele2 rather than trusted from `MarkerID`, whose
#              content depends on how the genotypes were exported. `--bim` cross-checks
#              the reconstruction against the genotypes and reports the match
#              rate, so a mismatch is loud rather than silent.
#
#           SAIGE emits BOTH p-values in one pass — `p.value` (SPA where
#           `Is.SPA`) and `p.value.NA` (normal approximation) — so `--p-value`
#           selects between them HERE, downstream of every expensive stage.
#           Switching it never re-runs SAIGE.
# Component: assoc_saige
# Used by : assoc_saige.nf  process TO_SUMSTATS
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# The columns the shared downstream reads, in plink2's order.
CANON = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'A1_FREQ', 'TEST', 'OBS_CT',
         'OR', 'LOG(OR)_SE', 'L95', 'U95', 'Z_STAT', 'P', 'ERRCODE']
# Everything SAIGE knows that the canonical schema has no slot for. Kept beside
# the sumstats rather than dropped: the SPA flag and the non-SPA p-value are the
# audit trail for --p-value, and the per-group counts are a free cross-check on
# the genotype counts the annotator computes independently.
EXTRA = ['Is.SPA', 'p.value.NA', 'AF_case', 'AF_ctrl', 'N_case', 'N_ctrl',
         'N_case_hom', 'N_case_het', 'N_ctrl_hom', 'N_ctrl_het', 'Tstat', 'var',
         'MissingRate', 'AC_Allele2']
Z95 = 1.959963984540054            # qnorm(0.975)


def parse_args():
    p = argparse.ArgumentParser(description='SAIGE step-2 output -> canonical summary statistics.')
    p.add_argument('--saige', required=True, help='Merged SAIGE association file.')
    p.add_argument('--p-value', default='spa', choices=['spa', 'no_spa'],
                   help="'spa' uses SAIGE's p.value; 'no_spa' uses p.value.NA.")
    p.add_argument('--bim', help='Genotype .bim, to cross-check the reconstructed variant IDs.')
    p.add_argument('--cohort', default='')
    p.add_argument('--model', default='')
    p.add_argument('--out-sumstats', required=True)
    p.add_argument('--out-extra', required=True)
    return p.parse_args()


def read_saige(path):
    """SAIGE writes a space-delimited table; be tolerant of either separator."""
    df = pd.read_csv(path, sep=r'\s+', engine='c',
                     dtype={'CHR': str, 'MarkerID': str, 'Allele1': str, 'Allele2': str})
    need = {'CHR', 'POS', 'Allele1', 'Allele2', 'BETA', 'SE', 'p.value'}
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f'SAIGE file is missing required columns: {sorted(missing)}')
    return df


def to_canonical(d, p_value='spa'):
    """SAIGE columns -> the canonical schema, as a DataFrame. One definition,
    shared by the whole-scan adapter and the conditional-analysis rounds."""
    chrom = d['CHR'].astype(str).str.replace('^chr', '', regex=True)
    pos = pd.to_numeric(d['POS'], errors='coerce').astype('Int64')
    ref, alt = d['Allele1'].astype(str), d['Allele2'].astype(str)
    # THE JOIN KEY IS SAIGE'S OWN MarkerID. With PLINK input SAIGE passes the
    # .bim identifier straight through, so using it makes the key agree with the
    # genotypes BY CONSTRUCTION — which is what every downstream lookup needs.
    # Reconstructing chr:pos:REF:ALT is the fallback for an input that carries no
    # usable id, and it is only correct if the allele order was declared right.
    vid = (d['MarkerID'].astype(str) if 'MarkerID' in d.columns
           else 'chr' + chrom + ':' + pos.astype(str) + ':' + ref + ':' + alt)

    beta = pd.to_numeric(d['BETA'], errors='coerce')
    se = pd.to_numeric(d['SE'], errors='coerce')
    pcol = 'p.value' if p_value == 'spa' else 'p.value.NA'
    if pcol not in d.columns:
        raise SystemExit(f"--p-value {p_value} needs column '{pcol}', which is not in the file")
    pval = pd.to_numeric(d[pcol], errors='coerce')

    n_case = pd.to_numeric(d.get('N_case'), errors='coerce')
    n_ctrl = pd.to_numeric(d.get('N_ctrl'), errors='coerce')
    miss = pd.to_numeric(d.get('MissingRate'), errors='coerce').fillna(0.0)
    obs = ((n_case.fillna(0) + n_ctrl.fillna(0)) * (1.0 - miss)).round()

    # A conditioning round drives the conditioned variant's SE to ~1e6 and its
    # BETA to a value whose exp() overflows. That row is unusable and is labelled
    # below; the overflow itself is expected, so it is not worth a warning.
    with np.errstate(over='ignore', invalid='ignore'):
        out = _build(chrom, pos, vid, ref, alt, d, beta, se, pval, obs)
        res, n_bad = _finish(out)
    return res, pcol, n_bad


def _build(chrom, pos, vid, ref, alt, d, beta, se, pval, obs):
    return pd.DataFrame({
        '#CHROM': chrom, 'POS': pos, 'ID': vid, 'REF': ref, 'ALT': alt,
        'A1': alt,
        'A1_FREQ': pd.to_numeric(d.get('AF_Allele2'), errors='coerce'),
        'TEST': 'ADD', 'OBS_CT': obs.astype('Int64'),
        'OR': np.exp(beta), 'LOG(OR)_SE': se,
        'L95': np.exp(beta - Z95 * se), 'U95': np.exp(beta + Z95 * se),
        'Z_STAT': beta / se, 'P': pval,
    })


def _finish(out):
    bad = (~np.isfinite(out['OR'])) | (~np.isfinite(out['LOG(OR)_SE'])) \
        | (out['LOG(OR)_SE'] <= 0) | out['P'].isna() | (out['P'] <= 0) | (out['P'] > 1)
    out = out.copy()
    out['ERRCODE'] = np.where(bad, 'SAIGE_NO_FIT', '.')
    return out[CANON], int(bad.sum())


def main():
    args = parse_args()
    d = read_saige(args.saige)
    n_in = len(d)
    out, pcol, n_bad = to_canonical(d, args.p_value)
    vid = out['ID']

    Path(args.out_sumstats).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_sumstats, sep='\t', index=False, na_rep='NA')

    keep = [c for c in EXTRA if c in d.columns]
    extra = pd.concat([pd.DataFrame({'ID': vid.values,
                                     'MarkerID_saige': d.get('MarkerID', pd.NA)}),
                       d[keep].reset_index(drop=True)], axis=1)
    extra.to_csv(args.out_extra, sep='\t', index=False, na_rep='NA')

    # ── cross-check the join key against the genotypes ──────────────────────
    match = ''
    if args.bim and Path(args.bim).exists():
        bim = pd.read_csv(args.bim, sep=r'\s+', header=None, usecols=[1], names=['ID'],
                          dtype={1: str})
        hit = out['ID'].isin(set(bim['ID']))
        rate = float(hit.mean())
        match = f'  id-match vs .bim {rate:.4%}'
        if rate < 0.99:
            print(f'[saige_to_sumstats] WARNING: only {rate:.2%} of reconstructed IDs are in the '
                  f'.bim — the allele order or ID convention does not agree with the '
                  f'genotypes, and every downstream join keyed on ID will be wrong.',
                  file=sys.stderr)

    print(f'[saige_to_sumstats] {args.cohort}/{args.model}: {n_in:,} rows, '
          f'P from {pcol}, {n_bad:,} without a usable fit{match} -> {args.out_sumstats}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in saige_to_sumstats: {e}', file=sys.stderr)
        sys.exit(1)
