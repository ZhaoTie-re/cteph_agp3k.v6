#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The credible-set variant table — one row per SuSiE credible-set
#           MEMBER, which is the unit downstream annotation work needs.
#           Until now the members existed only as a comma-joined string in
#           <peak>.cs.tsv and as a `cs` column in <peak>.pip.tsv; neither is a
#           usable input to anything else.
#
#           Each member carries its set-level context (size, coverage, purity),
#           its own posterior weight and rank within the set, and the three
#           annotations: rsID, nearest/overlapping gene, and snpEff functional
#           consequence from the platform's pre-computed index.
# Component: assoc_plink2
# Used by : assoc_plink2.nf  process CS_VARIANTS
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import variant_annot as VA

OUT_COLS = ['cohort', 'peak_id', 'cs', 'cs_size', 'cs_coverage', 'cs_purity_min_abs_corr',
            'variant_id', 'rsID', 'Gene', 'Gene_Biotype', 'chrom', 'pos', 'EA', 'OA',
            'A1_FREQ', 'P', 'OR', 'L95', 'U95', 'Beta', 'SE', 'pip', 'pip_rank', 'cum_pip',
            'is_lead', 'is_top_pip'] + VA.SNPEFF_FIELDS


def parse_args():
    p = argparse.ArgumentParser(description='Per-variant table of SuSiE credible-set members.')
    p.add_argument('--pip', required=True)
    p.add_argument('--cs', required=True)
    p.add_argument('--sumstat', help='<peak>.sumstat.tsv, for REF/ALT of each member.')
    p.add_argument('--cohort', required=True)
    p.add_argument('--peak-id', required=True)
    p.add_argument('--lead-id', required=True)
    p.add_argument('--tabix', required=True)
    p.add_argument('--rsid-vcf')
    p.add_argument('--snpeff-index')
    p.add_argument('--rscript')
    p.add_argument('--gene-script')
    p.add_argument('--out', required=True)
    return p.parse_args()


def main():
    args = parse_args()
    out_path = Path(args.out)
    work = out_path.parent / '_work'

    pip = pd.read_csv(args.pip, sep='\t', dtype={'ID': str, 'CHROM': str})
    cs = (pd.read_csv(args.cs, sep='\t')
          if Path(args.cs).exists() and Path(args.cs).stat().st_size else pd.DataFrame())

    members = pip[pip['cs'].notna()].copy() if 'cs' in pip else pip.iloc[0:0].copy()
    if not len(members) or not len(cs):
        pd.DataFrame(columns=OUT_COLS).to_csv(out_path, sep='\t', index=False)
        print(f'[cs_variants] {args.cohort} {args.peak_id}: no credible set')
        return

    # REF/ALT so the effect allele can be paired with the other allele, and so the
    # snpEff / rsID keys can be reconstructed exactly. The same file carries L95/U95,
    # which pip.tsv does not: an OR without its interval is not a reportable
    # estimate, so the CI is picked up here rather than left out of the table.
    ref_alt, ci = {}, {}
    if args.sumstat and Path(args.sumstat).exists():
        ss = pd.read_csv(args.sumstat, sep='\t', dtype={'ID': str})
        if {'REF', 'ALT'} <= set(ss.columns):
            ref_alt = {r['ID']: (r['REF'], r['ALT']) for _, r in ss.iterrows()}
        if {'L95', 'U95'} <= set(ss.columns):
            ci = {r['ID']: (r['L95'], r['U95']) for _, r in ss.iterrows()}

    meta = {int(r['cs']): r for _, r in cs.iterrows()}
    members['cs'] = members['cs'].astype(int)
    # Rank and cumulative mass are WITHIN a set, PIP-descending — the order any
    # downstream prioritisation will read.
    members = members.sort_values(['cs', 'pip'], ascending=[True, False])
    members['pip_rank'] = members.groupby('cs').cumcount() + 1
    members['cum_pip'] = members.groupby('cs')['pip'].cumsum()

    keys = [(r['ID'], r['CHROM'], int(r['POS'])) for _, r in members.iterrows()]
    rsids = VA.lookup_rsids(args.tabix, args.rsid_vcf, keys)
    snpeff = VA.lookup_snpeff(args.tabix, args.snpeff_index, keys)
    genes = VA.lookup_genes(args.rscript, args.gene_script,
                            [(r['CHROM'], int(r['POS'])) for _, r in members.iterrows()], work)

    rows = []
    for _, r in members.iterrows():
        vid = r['ID']
        m = meta.get(int(r['cs']), {})
        ea = str(r.get('A1', ''))
        ref, alt = ref_alt.get(vid, ('', ''))
        oa = (ref if ea == alt else alt) if (ref or alt) else ''
        orv = float(r['OR']) if pd.notna(r.get('OR')) else np.nan
        lo, hi = ci.get(vid, (np.nan, np.nan))
        gi = genes.get((str(r['CHROM']), int(r['POS'])))
        row = {
            'cohort': args.cohort, 'peak_id': args.peak_id, 'cs': int(r['cs']),
            'cs_size': int(m.get('size', 0)) or len(members[members.cs == r['cs']]),
            'cs_coverage': m.get('coverage', np.nan),
            'cs_purity_min_abs_corr': m.get('purity_min_abs_corr', np.nan),
            'variant_id': vid, 'rsID': rsids.get(vid, '.'),
            'Gene': (gi['Gene'] if gi is not None else '.'),
            'Gene_Biotype': (gi['Gene_Biotype'] if gi is not None else '.'),
            'chrom': r['CHROM'], 'pos': int(r['POS']), 'EA': ea, 'OA': oa,
            'A1_FREQ': r.get('A1_FREQ', np.nan), 'P': r.get('P', np.nan),
            'OR': round(orv, 5) if np.isfinite(orv) else np.nan,
            'L95': round(float(lo), 5) if pd.notna(lo) else np.nan,
            'U95': round(float(hi), 5) if pd.notna(hi) else np.nan,
            'Beta': round(float(np.log(orv)), 5) if np.isfinite(orv) and orv > 0 else np.nan,
            'SE': r.get('SE', np.nan),
            'pip': round(float(r['pip']), 6), 'pip_rank': int(r['pip_rank']),
            'cum_pip': round(float(r['cum_pip']), 6),
            'is_lead': bool(vid == args.lead_id),
            'is_top_pip': bool(int(r['pip_rank']) == 1),
        }
        row.update(snpeff.get(vid, VA.SNPEFF_MISSING))
        rows.append(row)

    df = pd.DataFrame(rows, columns=OUT_COLS)
    df.to_csv(out_path, sep='\t', index=False)
    for f in list(work.glob('annot_*.tsv')) if work.exists() else []:
        f.unlink()

    print(f'[cs_variants] {args.cohort} {args.peak_id}: {len(df)} member(s) in '
          f"{df['cs'].nunique()} credible set(s); rsID {int((df.rsID != '.').sum())}, "
          f"snpEff {int((df.snpEff_Effect != '.').sum())}, "
          f"gene {int((df.Gene != '.').sum())}")
    show = ['cs', 'pip_rank', 'variant_id', 'rsID', 'Gene', 'pip', 'cum_pip',
            'snpEff_Effect', 'snpEff_Impact']
    print(df[show].to_string(index=False))


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in cs_variants: {e}', file=sys.stderr)
        sys.exit(1)
