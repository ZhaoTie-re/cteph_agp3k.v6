#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Report EVERY additive peak lead in EVERY cohort, whether or not that
#           cohort called a peak there.
#
#           The three sample sets are nested (narrow subset of intermediate
#           subset of full), so a variant that is a peak in one of them always
#           has a perfectly good estimate in the other two. Reporting only the
#           cohort that happened to call the peak throws that away and leaves a
#           blank that is indistinguishable from a missing estimate — which is
#           exactly what the old comparison figure did (a lead called in two
#           rows, not three).
#
#           `called_peak` carries the distinction explicitly: genome_wide /
#           suggestive / not_a_peak. It is never inferred from a blank cell.
#
#           This is DESCRIPTION, not replication. See docs/METHODS.md §12.
# Component: assoc_plink2
# Used by : assoc_plink2.nf  process CROSS_COHORT
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import variant_annot as VA

USECOLS = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'A1_FREQ',
           'OBS_CT', 'OR', 'LOG(OR)_SE', 'L95', 'U95', 'P', 'ERRCODE']
OUT_COLS = ['variant_id', 'chrom', 'pos', 'best_tier',
            'cohort', 'called_peak', 'peak_id'] + VA.STAT_FIELDS
TIER_RANK = {'genome_wide': 0, 'suggestive': 1, 'not_a_peak': 2}


def parse_args():
    p = argparse.ArgumentParser(description='Every peak lead, in every cohort.')
    p.add_argument('--peaks', action='append', required=True, metavar='COHORT=PATH')
    p.add_argument('--glm', action='append', required=True, metavar='COHORT=PATH')
    p.add_argument('--bfile', action='append', required=True, metavar='COHORT=PREFIX')
    p.add_argument('--cohort-order', required=True)
    p.add_argument('--plink2', required=True)
    p.add_argument('--tabix', required=True)
    p.add_argument('--rsid-vcf')
    p.add_argument('--rscript')
    p.add_argument('--gene-script')
    p.add_argument('--snpeff-index')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--out', required=True)
    return p.parse_args()


def spec(items):
    """['a=1', 'b=2'] -> {'a': '1', 'b': '2'}"""
    d = {}
    for s in items:
        k, _, v = s.partition('=')
        if v:
            d[k] = v
    return d


def main():
    args = parse_args()
    out_path = Path(args.out)
    order = [c for c in args.cohort_order.split(',') if c]
    peaks_of, glm_of, bfile_of = spec(args.peaks), spec(args.glm), spec(args.bfile)
    order = [c for c in order if c in peaks_of and c in glm_of and c in bfile_of]

    # ── the union of every peak lead, and what each cohort made of it ────────
    per_cohort = {}
    for c in order:
        p = Path(peaks_of[c])
        d = (pd.read_csv(p, sep='\t', dtype={'chrom': str, 'lead_id': str})
             if p.exists() and p.stat().st_size else pd.DataFrame())
        per_cohort[c] = d
    called = {}                       # (cohort, variant_id) -> (tier, peak_id)
    union = {}                        # variant_id -> (chrom, pos)
    for c, d in per_cohort.items():
        for _, r in d.iterrows():
            vid = str(r['lead_id'])
            union.setdefault(vid, (str(r['chrom']), int(r['lead_pos'])))
            called[(c, vid)] = (r['tier'], r['peak_id'])
    if not union:
        pd.DataFrame(columns=OUT_COLS).to_csv(out_path, sep='\t', index=False)
        print('[cross_cohort] no peak leads in any cohort')
        return

    best_tier = {}
    for vid in union:
        tiers = [called[(c, vid)][0] for c in order if (c, vid) in called]
        best_tier[vid] = min(tiers, key=lambda t: TIER_RANK.get(t, 9)) if tiers else 'not_a_peak'

    # ── per cohort: pull those variants out of its own additive scan ─────────
    rows = []
    for c in order:
        head = pd.read_csv(glm_of[c], sep='\t', nrows=0)
        cols = [x for x in USECOLS if x in head.columns]
        glm = pd.read_csv(glm_of[c], sep='\t', usecols=cols,
                          dtype={'#CHROM': str, 'ID': str, 'ERRCODE': str})
        glm = glm.rename(columns={'#CHROM': 'CHROM', 'LOG(OR)_SE': 'SE'})
        glm = glm[glm['ID'].isin(union)]
        for x in ('POS', 'P', 'OR', 'SE', 'L95', 'U95', 'A1_FREQ', 'OBS_CT'):
            if x in glm:
                glm[x] = pd.to_numeric(glm[x], errors='coerce')

        recs = [{'variant_id': r['ID'], 'chrom': str(r['CHROM']), 'pos': int(r['POS']),
                 'a1': r.get('A1', ''), 'ref': r.get('REF', ''), 'alt': r.get('ALT', ''),
                 'or': r.get('OR'), 'se': r.get('SE'), 'l95': r.get('L95'),
                 'u95': r.get('U95'), 'p': r.get('P'), 'a1_freq': r.get('A1_FREQ'),
                 'obs_ct': r.get('OBS_CT')} for _, r in glm.iterrows()]
        annot = VA.annotate_variants(
            recs, bfile=bfile_of[c], plink2=args.plink2, tabix=args.tabix,
            work=out_path.parent / f'_work_{c}', threads=args.threads,
            rsid_vcf=args.rsid_vcf, rscript=args.rscript,
            gene_script=args.gene_script, snpeff_index=args.snpeff_index)

        for _, a in annot.iterrows():
            vid = a['variant_id']
            tier, pid = called.get((c, vid), ('not_a_peak', '.'))
            row = {'variant_id': vid, 'chrom': a['chrom'], 'pos': a['pos'],
                   'best_tier': best_tier.get(vid, 'not_a_peak'),
                   'cohort': c, 'called_peak': tier, 'peak_id': pid}
            row.update({k: a[k] for k in VA.STAT_FIELDS})
            rows.append(row)
        missing = set(union) - set(annot['variant_id']) if len(annot) else set(union)
        if missing:
            # A variant absent from a cohort's call set is a real fact about that
            # cohort, not a lookup failure — record it rather than dropping the row.
            for vid in sorted(missing):
                ch, po = union[vid]
                row = {'variant_id': vid, 'chrom': ch, 'pos': po,
                       'best_tier': best_tier.get(vid, 'not_a_peak'),
                       'cohort': c, 'called_peak': 'not_in_call_set', 'peak_id': '.'}
                row.update({k: ('.' if k in ('rsID', 'Gene', 'Gene_Biotype', 'EA', 'OA')
                                else np.nan) for k in VA.STAT_FIELDS})
                rows.append(row)
        print(f'[cross_cohort] {c}: {len(annot)} of {len(union)} lead variants present')

    # Genome-wide leads first, then by position, with a variant's three cohorts
    # adjacent and in nesting order — so the file reads the way it is used.
    out = pd.DataFrame(rows, columns=OUT_COLS)
    out['_t'] = out['best_tier'].map(TIER_RANK).fillna(9)
    out['_c'] = out['cohort'].map({c: i for i, c in enumerate(order)})
    out['_k'] = pd.to_numeric(out['chrom'].astype(str).str.replace('chr', '', regex=False),
                              errors='coerce').fillna(99)
    out = out.sort_values(['_t', '_k', 'pos', '_c']).drop(columns=['_t', '_c', '_k'])
    out.to_csv(out_path, sep='\t', index=False)

    n_gw = out[out.best_tier == 'genome_wide']['variant_id'].nunique()
    print(f'[cross_cohort] {out["variant_id"].nunique()} lead variants x {len(order)} cohorts '
          f'= {len(out)} rows ({n_gw} genome-wide somewhere)')
    if n_gw:
        show = ['variant_id', 'rsID', 'Gene', 'cohort', 'called_peak', 'OR', 'L95', 'U95',
                'P', 'Case_EAF', 'Control_EAF']
        print(out[out.best_tier == 'genome_wide'][show].to_string(index=False))


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in cross_cohort: {e}', file=sys.stderr)
        sys.exit(1)
