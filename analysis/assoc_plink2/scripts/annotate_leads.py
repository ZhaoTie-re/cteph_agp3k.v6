#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The two variant tables for one cohort.
#
#             lead_annotation.tsv        one row per ADDITIVE peak lead, both tiers
#             model_peaks_annotation.tsv  one row per peak of EVERY model, both tiers
#
#           Both carry the identical field set (variant_annot.STAT_FIELDS), so
#           they stack and a reader learns the columns once. The second table is
#           what each scan figure draws itself from: before it existed, DOM and
#           REC peaks were a bare count in scan_qc.tsv and their Manhattans
#           carried ADD's labels or none.
#
#           Every variant is annotated in ONE pass over the union of the two
#           sets — the two plink2 --geno-counts/--hardy calls just get a longer
#           --extract list, so a second table costs nothing.
# Component: assoc_plink2
# Used by : assoc_plink2.nf  process ANNOTATE_LEADS
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import variant_annot as VA

LEAD_COLS = ['cohort', 'peak_id', 'tier', 'chrom', 'pos', 'variant_id'] + VA.STAT_FIELDS
HIT_COLS = ['cohort', 'model', 'peak_id', 'tier', 'n_sig_variants',
            'n_genomewide_variants', 'chrom', 'pos', 'variant_id'] + VA.STAT_FIELDS


def parse_args():
    p = argparse.ArgumentParser(description='Annotate this cohort\'s peak leads and genome-wide hits.')
    p.add_argument('--peaks', required=True)
    p.add_argument('--model-peaks', help='model_peaks.tsv from call_peaks — all three models.')
    p.add_argument('--bfile', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--plink2', required=True)
    p.add_argument('--tabix', required=True)
    p.add_argument('--rsid-vcf', help='Tabix-indexed VCF carrying rsIDs; omit to leave rsID as "."')
    p.add_argument('--rscript', help='Rscript in the R env; omit to leave Gene as "."')
    p.add_argument('--gene-script', help='Path to gene_annotate.R')
    p.add_argument('--snpeff-index', help='Pre-computed snpEff index directory')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--out', required=True)
    p.add_argument('--out-model-peaks', required=True)
    return p.parse_args()


def record(r, prefix='lead_'):
    """One association record in the shape variant_annot.annotate_variants wants."""
    return {'variant_id': r[f'{prefix}id'], 'chrom': str(r['chrom']),
            'pos': int(r[f'{prefix}pos'] if f'{prefix}pos' in r else r['pos']),
            'a1': r.get(f'{prefix}a1', ''), 'ref': r.get(f'{prefix}ref', ''),
            'alt': r.get(f'{prefix}alt', ''), 'or': r.get(f'{prefix}or'),
            'se': r.get(f'{prefix}se'), 'l95': r.get(f'{prefix}l95'),
            'u95': r.get(f'{prefix}u95'), 'p': r.get(f'{prefix}p'),
            'a1_freq': r.get(f'{prefix}a1_freq'), 'obs_ct': r.get(f'{prefix}obs_ct')}


def read(path, **kw):
    if not path or not Path(path).exists() or not Path(path).stat().st_size:
        return pd.DataFrame()
    return pd.read_csv(path, sep='\t', **kw)


def main():
    args = parse_args()
    out_path, hits_path = Path(args.out), Path(args.out_model_peaks)
    work = out_path.parent / '_work'

    peaks = read(args.peaks, dtype={'chrom': str, 'lead_id': str})
    hits = read(args.model_peaks, dtype={'chrom': str, 'lead_id': str})

    peak_recs = [record(r) for _, r in peaks.iterrows()] if len(peaks) else []
    hit_recs = [record(r) for _, r in hits.iterrows()] if len(hits) else []
    if not peak_recs and not hit_recs:
        pd.DataFrame(columns=LEAD_COLS).to_csv(out_path, sep='\t', index=False)
        pd.DataFrame(columns=HIT_COLS).to_csv(hits_path, sep='\t', index=False)
        print(f'[annotate_leads] {args.cohort}: no peaks in any model')
        return

    # The two tables are annotated SEPARATELY, on purpose. A variant can be a
    # genome-wide hit under one model and a peak lead under another — chr5:66174098
    # is a dominant hit at P = 4.6e-8 and an additive suggestive lead at P = 1.7e-7 —
    # so a single union keyed by variant_id would give one of them the other's
    # effect and P. Reporting a model's result under another model's statistics is
    # exactly what METHODS §4 rules out. The variant-level lookups repeat for the
    # overlap, which costs one extra pair of plink2 reports per cohort.
    def assemble(src, recs, prefix_cols, cols, tag):
        if not recs:
            return pd.DataFrame(columns=cols)
        annot = VA.annotate_variants(
            recs, bfile=args.bfile, plink2=args.plink2, tabix=args.tabix,
            work=work / tag, threads=args.threads, rsid_vcf=args.rsid_vcf,
            rscript=args.rscript, gene_script=args.gene_script,
            snpeff_index=args.snpeff_index)
        # Zip POSITIONALLY, never by variant_id. In model_peaks the same variant
        # appears once per model — a locus can be a genome-wide peak under more
        # dominant suggestive one — so a {variant_id: row} map keeps whichever
        # model was annotated last and hands its OR to the others. That silently
        # gave a peak of one model the odds ratio fitted under another.
        # annotate_variants appends one output row per input record, in order, so
        # position is the only key that distinguishes them.
        assert len(annot) == len(recs), f'{tag}: {len(annot)} annotated vs {len(recs)} records'
        rows = []
        for (_, r), (_, a) in zip(src.iterrows(), annot.iterrows()):
            row = {k: r[k] for k in prefix_cols}
            row.update({k: a[k] for k in ['chrom', 'pos', 'variant_id'] + VA.STAT_FIELDS})
            rows.append(row)
        return pd.DataFrame(rows, columns=cols)

    leads = assemble(peaks, peak_recs, ['cohort', 'peak_id', 'tier'], LEAD_COLS, 'leads')
    gwh = assemble(hits, hit_recs,
                   ['cohort', 'model', 'peak_id', 'tier', 'n_sig_variants',
                    'n_genomewide_variants'], HIT_COLS, 'model_peaks')

    # A check that can fail: every row's P must match the tier it claims. This is
    # what catches an association statistic that belongs to a different scan —
    # chr5:66174098 is a dominant genome-wide peak AND an additive suggestive one,
    # and an earlier union-by-variant gave the dominant row additive statistics.
    for frame, name in ((gwh, 'model_peaks'), (leads, 'lead')):
        if not len(frame):
            continue
        bad = frame[((frame.tier == 'genome_wide') & (frame.P >= 5e-8))
                    | ((frame.tier == 'suggestive') & ((frame.P >= 1e-5) | (frame.P < 5e-8)))]
        if len(bad):
            print(f'ERROR: {len(bad)} {name} row(s) carry a P that contradicts their tier:\n'
                  f'{bad[["tier", "P"]].head(10)}', file=sys.stderr)
            sys.exit(1)
    leads.to_csv(out_path, sep='\t', index=False)
    gwh.to_csv(hits_path, sep='\t', index=False)

    n_gw = int((leads.tier == 'genome_wide').sum()) if len(leads) else 0
    print(f'[annotate_leads] {args.cohort}: {len(leads)} additive peak leads '
          f'({n_gw} genome-wide, {len(leads) - n_gw} suggestive), '
          f'{len(gwh)} peak(s) across all models; '
          f"rsID {int((leads.rsID != '.').sum()) if len(leads) else 0}, "
          f"gene {int((leads.Gene != '.').sum()) if len(leads) else 0}, "
          f"snpEff {int((leads.snpEff_Effect != '.').sum()) if len(leads) else 0}")
    show = ['peak_id', 'tier', 'rsID', 'Gene', 'EA', 'OA', 'OR', 'L95', 'U95', 'P',
            'Case_Genotype_Distribution', 'Case_EAF', 'Control_EAF']
    if len(leads):
        print(leads[leads.tier == 'genome_wide'][show].to_string(index=False))
    if len(gwh):
        g = gwh[gwh.tier == 'genome_wide']
        if len(g):
            print(g[['model', 'peak_id', 'rsID', 'Gene', 'OR', 'L95', 'U95', 'P']]
                  .to_string(index=False))
        print('   peaks per model/tier:',
              gwh.groupby(['model', 'tier']).size().to_dict())


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in annotate_leads: {e}', file=sys.stderr)
        sys.exit(1)
