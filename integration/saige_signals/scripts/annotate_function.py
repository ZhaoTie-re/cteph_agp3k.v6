#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Functional annotation of every table variant.
#             snpEff  the JHRP v6 pre-computed index (as analysis/assoc_gene),
#                     most severe transcript + every consequence term + LoF/NMD
#             gene    Ensembl 86 overlapping / nearest gene with signed distance
#                     (gene_annotate.R, staged beside gene_utils.R)
# Inputs  : variants.tsv
# Outputs : function.tsv  one row per variant
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sig_common import chrom_plain, import_annot, load_variants, log, run, write_tsv  # noqa: E402


def fetch_snpeff(v, tabix, index_dir, va):
    keys = [(r.variant_id, r.chrom, r.pos) for _, r in v.iterrows()]
    base = va.lookup_snpeff(tabix, str(index_dir), keys)
    if len(base) != len(v):
        log(f'    snpEff: no index row for {len(v) - len(base)} variant(s): '
            + ', '.join(sorted(set(v.variant_id) - set(base))))
    empty = {k: '.' for k in ('snpEff_Effect', 'snpEff_Impact', 'snpEff_Gene', 'snpEff_FeatureID',
                              'snpEff_Biotype', 'snpEff_HGVS_c', 'snpEff_HGVS_p', 'snpEff_all_effects')}
    empty['snpEff_n_annotations'] = 0
    # LoF / NMD flags live in columns 26-33 of the index, which lookup_snpeff
    # does not carry; read them from the same rows.
    lof = {}
    for vid, chrom, pos in keys:
        path = Path(index_dir) / va.SNPEFF_TEMPLATE.format(chrom=chrom_plain(chrom))
        res = run([tabix, str(path), f'chr{chrom_plain(chrom)}:{pos}-{pos}'])
        flags = []
        for line in res.splitlines():
            f = line.split('\t')
            if len(f) < 34 or f'chr{chrom_plain(f[0])}:{f[1]}:{f[2]}:{f[3]}' != vid:
                continue
            if f[26] not in ('.', ''):
                flags.append(f'LoF {f[26]} ({f[29]})')
            if f[30] not in ('.', ''):
                flags.append(f'NMD {f[30]} ({f[33]})')
        lof[vid] = '; '.join(sorted(set(flags))) if flags else ''
    rows = []
    for vid in v.variant_id:
        b = base.get(vid, empty)
        rows.append({'variant_id': vid, 'snpEff_Effect': b['snpEff_Effect'], 'snpEff_Impact': b['snpEff_Impact'],
                     'snpEff_Gene': b['snpEff_Gene'], 'snpEff_FeatureID': b['snpEff_FeatureID'],
                     'snpEff_Biotype': b['snpEff_Biotype'], 'snpEff_HGVS_c': b['snpEff_HGVS_c'],
                     'snpEff_HGVS_p': b['snpEff_HGVS_p'], 'snpEff_n_annotations': int(b['snpEff_n_annotations']),
                     'snpEff_all_effects': b['snpEff_all_effects'], 'snpEff_LoF_NMD': lof[vid]})
    return pd.DataFrame(rows), len(base)


def fetch_genes(v, rscript, gene_script, va):
    """Ensembl 86 nearest/overlapping gene with signed distance, for all rows
    (cs_variants_all carries the symbol but not the distance)."""
    genes = va.lookup_genes(rscript, gene_script, [(r.chrom, r.pos) for _, r in v.iterrows()], Path('genes'))
    # lookup_genes returns {} silently when R fails or gene_utils.R is not beside gene_annotate.R
    assert len(genes) == len(v), f'gene annotation for {len(genes)}/{len(v)} positions'
    rows = []
    for _, r in v.iterrows():
        g = genes[(str(r.chrom), int(r.pos))]
        d = str(g['Gene_Distance_bp'])
        rows.append({'variant_id': r.variant_id, 'Gene_EnsDb': g['Gene'], 'Gene_Biotype': g['Gene_Biotype'],
                     'Gene_Distance_bp': int(float(d)) if d not in ('.', '', 'nan') else 0})
    out = pd.DataFrame(rows)
    chk = out.merge(v[['variant_id', 'Gene']], on='variant_id')
    bad = chk[chk.Gene_EnsDb != chk.Gene]
    assert bad.empty, bad[['variant_id', 'Gene', 'Gene_EnsDb']].to_string()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--variants', required=True)
    ap.add_argument('--tabix', required=True)
    ap.add_argument('--snpeff-index', required=True, help='directory of all.VQSR3.chr<N>.vcf_out.tsv.2.gz')
    ap.add_argument('--rscript', required=True)
    ap.add_argument('--gene-script', required=True, help='gene_annotate.R, staged next to gene_utils.R')
    ap.add_argument('--annot-lib', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    assert Path(a.snpeff_index).is_dir(), a.snpeff_index     # lookup_snpeff returns {} silently otherwise
    va = import_annot(a.annot_lib)
    v = load_variants(a.variants)
    snpeff, n_hit = fetch_snpeff(v, a.tabix, a.snpeff_index, va)
    genes = fetch_genes(v, a.rscript, a.gene_script, va)
    out = snpeff.merge(genes, on='variant_id')
    assert len(out) == len(v)
    write_tsv(out, a.out)
    log(f'annotation: snpEff {n_hit}/{len(v)}, Ensembl gene {len(genes)}/{len(v)}; '
        f'impact {out.snpEff_Impact.value_counts().to_dict()}')


if __name__ == '__main__':
    main()
