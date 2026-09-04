#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Both engines -> the canonical 16-column plink2 --glm schema, so the
#           two models are read by exactly the same downstream code and the
#           comparison between them is a comparison of ESTIMATES, not of formats.
#
#           THE EFFECT ALLELE IS 'P', ALWAYS. build_hla_markers.py writes
#           REF = A (the HLA allele or residue is absent) and ALT = P (present),
#           verified through plink2 --freq: ALT_FREQS reproduces the P frequency
#           exactly. So a positive log-odds means CARRYING the allele raises
#           risk, at every marker, in both models. This is asserted here rather
#           than trusted, because an inverted effect allele silently inverts
#           every odds ratio in the paper.
#
#           assoc_saige has two drifts this component deliberately does not
#           inherit: --AlleleOrder is 'alt-first' in its code and 'ref-first' in
#           its METHODS, and its ID is MarkerID-first while documented as
#           reconstructed. Here the ID is the marker id we minted, both engines
#           are given the same .bim, and the ID is checked against it.
# Component: assoc_hla
# ---------------------------------------------------------------------------
import argparse
import sys

import numpy as np
import pandas as pd

CANON = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'A1_FREQ', 'TEST', 'OBS_CT',
         'OR', 'LOG(OR)_SE', 'L95', 'U95', 'Z_STAT', 'P', 'ERRCODE']
Z95 = 1.959963984540054
SAIGE_EXTRA = ['Is.SPA', 'p.value.NA', 'AF_case', 'AF_ctrl', 'N_case', 'N_ctrl',
               'N_case_hom', 'N_case_het', 'N_ctrl_hom', 'N_ctrl_het', 'Tstat',
               'var', 'MissingRate', 'AC_Allele2']


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assoc', required=True, help='plink2 --glm or SAIGE step-2 output')
    p.add_argument('--engine', required=True, choices=('plink2', 'saige'))
    p.add_argument('--bim', required=True, help='the .bim both engines were given')
    p.add_argument('--marker-map', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--model', required=True, help='fixed | random')
    p.add_argument('--marker-class', required=True, choices=('allele', 'residue'))
    p.add_argument('--p-value', default='spa', choices=('spa', 'no_spa'),
                   help='saige only: p.value (SPA) or p.value.NA (normal approximation)')
    p.add_argument('--out-sumstats', default='sumstats.tsv')
    p.add_argument('--out-extra', default='engine_extra.tsv')
    return p.parse_args()


def _finish(out, bim, mm, args):
    with np.errstate(over='ignore', invalid='ignore'):
        bad = ((~np.isfinite(out['OR'])) | (~np.isfinite(out['LOG(OR)_SE']))
               | (out['LOG(OR)_SE'] <= 0) | out['P'].isna()
               | (out['P'] <= 0) | (out['P'] > 1))
    out['ERRCODE'] = np.where(bad, f'{args.engine.upper()}_NO_FIT', '.')

    # The effect allele must be P at every marker, or the odds ratios are inverted.
    wrong = out.loc[out.A1 != 'P', 'ID']
    if len(wrong):
        raise SystemExit(f'ABORT: {len(wrong)} marker(s) report A1 != "P", so their effect '
                         f'is per ABSENCE of the HLA allele. First: {list(wrong[:3])}. '
                         f'build_hla_markers.py writes REF=A/ALT=P; something re-derived '
                         f'the allele order.')
    ids = set(bim.ID)
    miss = [i for i in out.ID if i not in ids]
    if miss:
        raise SystemExit(f'ABORT: {len(miss)} result id(s) are not in the .bim, so every '
                         f'downstream join keyed on ID is wrong. First: {miss[:3]}')

    out = out.merge(mm[['id', 'gene', 'marker_class', 'position']],
                    left_on='ID', right_on='id', how='left').drop(columns='id')
    out['cohort'], out['model'] = args.cohort, args.model
    return out


def from_plink2(d, args):
    d = d.rename(columns={c: c.strip() for c in d.columns})
    need = {'#CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'OR', 'LOG(OR)_SE', 'P'}
    if not need.issubset(d.columns):
        raise SystemExit(f'ABORT: {args.assoc} lacks {sorted(need - set(d.columns))}')
    d = d[d.get('TEST', 'ADD') == 'ADD'] if 'TEST' in d.columns else d
    out = pd.DataFrame({c: d[c] if c in d.columns else pd.NA for c in CANON})
    out['#CHROM'] = out['#CHROM'].astype(str).str.replace('^chr', '', regex=True)
    out['TEST'] = 'ADD'
    for c in ('OR', 'LOG(OR)_SE', 'L95', 'U95', 'Z_STAT', 'P', 'A1_FREQ'):
        out[c] = pd.to_numeric(out[c], errors='coerce')
    return out, pd.DataFrame({'ID': d['ID']})


def from_saige(d, args):
    need = {'CHR', 'POS', 'MarkerID', 'Allele1', 'Allele2', 'BETA', 'SE'}
    if not need.issubset(d.columns):
        raise SystemExit(f'ABORT: {args.assoc} lacks {sorted(need - set(d.columns))}')
    pcol = 'p.value' if args.p_value == 'spa' else 'p.value.NA'
    if pcol not in d.columns:
        raise SystemExit(f'ABORT: --p-value {args.p_value} needs column {pcol!r}; '
                         f'present: {[c for c in d.columns if "p.value" in c]}')
    beta = pd.to_numeric(d['BETA'], errors='coerce')
    se = pd.to_numeric(d['SE'], errors='coerce')
    # SAIGE's BETA is per Allele2, and Allele2 is the ALT allele. plink2 --make-bed
    # from our VCF writes the .bim with ALT in column 5 -- confirmed by --freq,
    # where ALT_FREQS reproduces the P frequency exactly -- so step 2 is run with
    # --AlleleOrder=alt-first, which tells SAIGE column 5 is ALT. Allele2 is then
    # P and BETA is per CARRYING the allele. The A1 == 'P' assertion in _finish()
    # is what proves this run by run rather than leaving it to this comment.
    with np.errstate(over='ignore', invalid='ignore'):
        out = pd.DataFrame({
            '#CHROM': d['CHR'].astype(str).str.replace('^chr', '', regex=True),
            'POS': pd.to_numeric(d['POS'], errors='coerce').astype('Int64'),
            'ID': d['MarkerID'].astype(str),
            'REF': d['Allele1'], 'ALT': d['Allele2'], 'A1': d['Allele2'],
            'A1_FREQ': pd.to_numeric(d.get('AF_Allele2'), errors='coerce'),
            'TEST': 'ADD',
            'OBS_CT': ((pd.to_numeric(d.get('N_case'), errors='coerce').fillna(0)
                        + pd.to_numeric(d.get('N_ctrl'), errors='coerce').fillna(0))
                       * (1.0 - pd.to_numeric(d.get('MissingRate'), errors='coerce')
                          .fillna(0))).round().astype('Int64'),
            'OR': np.exp(beta), 'LOG(OR)_SE': se,
            'L95': np.exp(beta - Z95 * se), 'U95': np.exp(beta + Z95 * se),
            'Z_STAT': beta / se,
            'P': pd.to_numeric(d[pcol], errors='coerce'),
            'ERRCODE': '.'})
    extra = pd.DataFrame({'ID': d['MarkerID'].astype(str)})
    for c in SAIGE_EXTRA:
        if c in d.columns:
            extra[c] = d[c]
    return out, extra


def main():
    args = parse_args()
    d = pd.read_csv(args.assoc, sep=r'\s+', dtype=str, comment=None)
    bim = pd.read_csv(args.bim, sep=r'\s+', header=None,
                      names=['CHROM', 'ID', 'CM', 'POS', 'A1', 'A2'], dtype=str)
    mm = pd.read_csv(args.marker_map, sep='\t')

    out, extra = (from_plink2 if args.engine == 'plink2' else from_saige)(d, args)
    out = _finish(out[CANON].copy(), bim, mm, args)

    keep = out.marker_class == args.marker_class
    if not keep.any():
        raise SystemExit(f'ABORT: no {args.marker_class} marker in {args.assoc}')
    out = out[keep]
    out[CANON + ['gene', 'position', 'cohort', 'model']].to_csv(
        args.out_sumstats, sep='\t', index=False, na_rep='NA')
    extra[extra.ID.isin(set(out.ID))].to_csv(args.out_extra, sep='\t',
                                             index=False, na_rep='NA')

    n_bad = int((out.ERRCODE != '.').sum())
    print(f'[hla_to_sumstats] {args.cohort}/{args.model}/{args.marker_class} '
          f'[{args.engine}]: {len(out):,} marker(s) -> {args.out_sumstats}')
    print(f'    effect allele P at every marker (asserted); {n_bad} unusable fit(s)')
    ok = out[out.ERRCODE == '.']
    if len(ok):
        print(f'    min P {ok.P.min():.3e} at {ok.loc[ok.P.idxmin(), "ID"]}; '
              f'OR range {ok.OR.min():.3f}-{ok.OR.max():.3f}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in hla_to_sumstats: {e}', file=sys.stderr)
        sys.exit(1)
