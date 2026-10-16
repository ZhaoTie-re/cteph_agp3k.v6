#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : ANNOVAR filter-based annotations of every table variant, in one run
#           (the resources of the common-variant table):
#             gnomAD v4.1 genomes   AF, AF_eas, AF_nfe of the ALT allele
#             dbSNP 151             allele-specific rsID (avsnp151), used to attribute
#                                   the rsID at multi-allelic ToMMo records
#           ANNOVAR re-encodes indels; its coordinates are asserted to be what
#           our REF/ALT imply, so every row is the variant asked for.
# Inputs  : variants.tsv; the ANNOVAR directory and the protocol names
# Outputs : annovar.tsv  one row per variant
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gs_common import load_variants, log, num, read_tsv, run, write_tsv  # noqa: E402


def blank(x):
    return '' if x in ('.', '', None) else str(x)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--variants', required=True)
    ap.add_argument('--annovar-dir', required=True)
    ap.add_argument('--perl', required=True)
    ap.add_argument('--gnomad', default='gnomad41_genome')
    ap.add_argument('--dbsnp', default='avsnp151')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    annovar = Path(a.annovar_dir)
    for proto in (a.gnomad, a.dbsnp):
        db = annovar / 'humandb' / f'hg38_{proto}.txt'
        assert db.is_file(), db
        with open(db) as fh:                       # readable: ANNOVAR reports a permission error only in its log
            fh.read(1)
    v = load_variants(a.variants).drop_duplicates('variant_id').reset_index(drop=True)
    w = Path('annovar')
    w.mkdir(exist_ok=True)
    vcf = w / 'variants.vcf'
    with open(vcf, 'w') as fh:
        fh.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for _, r in v.sort_values(['chrom', 'pos'], key=lambda s: s.astype(int) if s.name != 'chrom'
                                  else s.astype(int)).iterrows():
            fh.write(f'chr{r.chrom}\t{r.pos}\t{r.variant_id}\t{r.OA}\t{r.EA}\t.\t.\t.\n')
    av = w / 'variants.avinput'
    av.write_text(run([a.perl, annovar / 'convert2annovar.pl', '-format', 'vcf4', '-includeinfo', vcf]))
    assert sum(1 for _ in open(av)) == len(v)
    protos = [a.gnomad, a.dbsnp]
    run([a.perl, annovar / 'table_annovar.pl', av, str(annovar / 'humandb') + '/', '-buildver', 'hg38',
         '-protocol', ','.join(protos), '-operation', ','.join('f' * len(protos)), '-nastring', '.',
         '-otherinfo', '-remove', '-out', w / 'gs'])
    ma = read_tsv(w / 'gs.hg38_multianno.txt')
    assert len(ma) == len(v), len(ma)
    other = [c for c in ma.columns if c.startswith('Otherinfo')]
    ids = set(v.variant_id)

    def vid_of(row):
        hits = [row[c] for c in other if row[c] in ids]
        assert len(hits) == 1, hits
        return hits[0]
    ma['variant_id'] = ma.apply(vid_of, axis=1)
    assert ma.variant_id.is_unique
    for _, r in ma.merge(v[['variant_id', 'pos', 'OA', 'EA']], on='variant_id').iterrows():
        ref, alt, pos = r.OA, r.EA, int(r.pos)
        if len(ref) == len(alt) == 1:
            exp = (pos, pos, ref, alt)
        elif len(ref) > len(alt) and ref.startswith(alt):        # deletion
            exp = (pos + len(alt), pos + len(ref) - 1, ref[len(alt):], '-')
        elif len(alt) > len(ref) and alt.startswith(ref):        # insertion
            exp = (pos + len(ref) - 1, pos + len(ref) - 1, '-', alt[len(ref):])
        else:
            exp = None
        got = (int(r.Start), int(r.End), r.Ref, r.Alt)
        assert exp is None or got == exp, (r.variant_id, got, exp)
    g = a.gnomad
    out = pd.DataFrame({
        'variant_id': ma.variant_id,
        'gnomad_af': num(ma[f'{g}_AF']), 'gnomad_eas_af': num(ma[f'{g}_AF_eas']),
        'gnomad_nfe_af': num(ma[f'{g}_AF_nfe']),
        'dbsnp_rsid': ma[a.dbsnp].map(blank)})
    write_tsv(out, a.out)
    log(f'ANNOVAR ({", ".join(protos)}): gnomAD {out.gnomad_af.notna().sum()}/{len(out)} '
        f'(EAS {out.gnomad_eas_af.notna().sum()}), dbSNP rsID {int((out.dbsnp_rsid != "").sum())}')


if __name__ == '__main__':
    main()
