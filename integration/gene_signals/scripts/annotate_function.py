#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Functional annotation of every table variant, in the columns of the
#           common-variant table (saige_signals), for the gene the variant was
#           tested in:
#             consequence  the tested gene's most severe snpEff term, chosen as the
#                          gene map chose it (lowest impact rank, then the term)
#                          and asserted equal to the map's
#             transcript   the tested gene's most severe transcript row (protein-
#                          coding first): transcript ID, biotype, HGVS c. / p.
#             other        every other consequence term at the variant, over all
#                          transcripts and features of all genes (as saige_signals)
#             LoF / NMD    snpEff's predictions over the variant's rows
#             gene         Ensembl 86 gene overlapping the variant, else the nearest
#                          (gene_annotate.R, staged beside gene_utils.R; the function
#                          behind the common-variant table)
#           A 'sequence_feature' term is a protein-feature annotation over the
#           position, not a coding change; the transcript row places the variant.
# Inputs  : variants.tsv
# Outputs : function.tsv  one row per gene x variant
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gs_common import import_annot, load_variants, log, run, write_tsv  # noqa: E402

TEMPLATE = 'all.VQSR3.chr{chrom}.vcf_out.tsv.2.gz'
IMPACT_RANK = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}
COL = {'chrom': 0, 'pos': 1, 'ref': 2, 'alt': 3, 'effect': 8, 'impact': 9, 'gene': 10, 'geneid': 11, 'feature': 12,
       'featureid': 13, 'biotype': 14, 'hgvs_c': 16, 'hgvs_p': 17, 'lof_gene': 26, 'lof_perc': 29, 'nmd_gene': 30,
       'nmd_perc': 33}


def rows_of(tabix, index_dir, chrom, pos, vid):
    path = Path(index_dir) / TEMPLATE.format(chrom=chrom)
    out = []
    for line in run([tabix, path, f'chr{chrom}:{pos}-{pos}']).splitlines():
        f = line.split('\t')
        if len(f) < 34:
            continue
        g = {k: f[i] for k, i in COL.items()}
        if f'chr{g["chrom"].replace("chr", "")}:{g["pos"]}:{g["ref"]}:{g["alt"]}' == vid:
            out.append(g)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--variants', required=True)
    ap.add_argument('--tabix', required=True)
    ap.add_argument('--snpeff-index', required=True, help='directory of all.VQSR3.chr<N>.vcf_out.tsv.2.gz')
    ap.add_argument('--rscript', required=True)
    ap.add_argument('--gene-script', required=True, help='gene_annotate.R, staged next to gene_utils.R')
    ap.add_argument('--annot-lib', required=True, help='the staged _shared/scripts/variant_annot.py')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    assert Path(a.snpeff_index).is_dir(), a.snpeff_index
    va = import_annot(a.annot_lib)
    v = load_variants(a.variants)
    out, n_tx = [], 0
    for _, r in v.iterrows():
        rows = rows_of(a.tabix, a.snpeff_index, r.chrom, r.pos, r.variant_id)
        mine = [g for g in rows if g['gene'] == r.gene_symbol]
        assert mine, f'{r.variant_id}: no snpEff row for {r.gene_symbol}'
        best = min(mine, key=lambda g: (IMPACT_RANK.get(g['impact'], 99), g['effect']))
        # the consequence the gene map assigned, reproduced
        assert (best['effect'], best['impact']) == (r.effect, r.impact), \
            (r.variant_id, best['effect'], best['impact'], r.effect, r.impact)
        tx = [g for g in mine if g['feature'] == 'transcript']
        txb = min(tx, key=lambda g: (IMPACT_RANK.get(g['impact'], 99), g['biotype'] != 'protein_coding',
                                     g['effect'], g['featureid'])) if tx else None
        n_tx += bool(txb)
        shown = set(best['effect'].split('&'))
        others = list(dict.fromkeys(t for g in sorted(rows, key=lambda g: IMPACT_RANK.get(g['impact'], 99))
                                    for t in g['effect'].split('&') if t not in shown and t not in ('', '.')))
        lof = sorted({f'LoF {g["lof_gene"]} ({g["lof_perc"]})' for g in rows if g['lof_gene'] not in ('', '.')}
                     | {f'NMD {g["nmd_gene"]} ({g["nmd_perc"]})' for g in rows if g['nmd_gene'] not in ('', '.')})
        out.append({'set_name': r.set_name, 'variant_id': r.variant_id, 'consequence': best['effect'],
                    'impact': best['impact'],
                    'tx_gene': txb['gene'] if txb else r.gene_symbol,
                    'transcript': txb['featureid'] if txb else '',
                    'transcript_consequence': txb['effect'] if txb else '',
                    'transcript_biotype': txb['biotype'] if txb else '',
                    'hgvs_c': txb['hgvs_c'] if txb and txb['hgvs_c'] not in ('', '.') else '',
                    'hgvs_p': txb['hgvs_p'] if txb and txb['hgvs_p'] not in ('', '.') else '',
                    'other_consequences': ','.join(others), 'lof_nmd': '; '.join(lof),
                    'n_annotations': len(rows),
                    'other_genes': ','.join(sorted({g['gene'] for g in rows if g['gene'] not in (r.gene_symbol, '.', '')}))})
    f = pd.DataFrame(out)
    assert len(f) == len(v)
    # Ensembl 86 gene, as in the common-variant table
    genes = va.lookup_genes(a.rscript, a.gene_script, [(str(r.chrom), int(r.pos)) for _, r in v.iterrows()],
                            Path('genes'))
    # lookup_genes returns {} silently when R fails or gene_utils.R is not beside gene_annotate.R
    assert len(genes) == v[['chrom', 'pos']].drop_duplicates().shape[0], f'Ensembl gene for {len(genes)} positions'
    ens = []
    for _, r in v.iterrows():
        g = genes[(str(r.chrom), int(r.pos))]
        d = str(g['Gene_Distance_bp'])
        ens.append({'Gene_EnsDb': g['Gene'], 'Gene_Biotype': g['Gene_Biotype'],
                    'Gene_Distance_bp': int(float(d)) if d not in ('.', '', 'nan') else 0})
    f = pd.concat([f, pd.DataFrame(ens)], axis=1)
    write_tsv(f, a.out)
    n_same = int((f.Gene_EnsDb == v.gene_symbol.values).sum())
    log(f'annotation: {len(f)} gene x variant rows, every consequence equal to the gene map; transcript row for '
        f'{n_tx}; impact {f.impact.value_counts().to_dict()}; LoF/NMD flagged {int((f.lof_nmd != "").sum())}; '
        f'Ensembl 86 gene = tested gene for {n_same}/{len(f)}')


if __name__ == '__main__':
    main()
