#!/usr/bin/env python3
"""One biallelic record per fixed site, for every sample of a joint call.

GATK 4.1.4.0 GenotypeGVCFs has no --force-output-intervals, so the joint call is made twice:
  --variants  variant sites over the calling intervals (-stand-call-conf 0: a site carried weakly by
              these samples alone is not lost to QUAL), split, left-aligned, IDs set;
  --allsites  --include-non-variant-sites at the fixed positions only: what a sample's reference
              block says where no sample of the call carries an alternate allele.

Per fixed ID (chrom:pos:ref:alt):
  1. the variant record with that ID, unchanged;
  2. else a record at the same position whose alternate allele is another one (these samples carry
     only that one): every sample 0/0 for the fixed allele, AD = (ref depth, 0), DP and GQ kept —
     what `bcftools norm -m -any` gives a sample of another allele;
  3. else the non-variant record: 0/0 with GQ <- RGQ and AD <- (DP, 0). Without AD the cohort's
     AlleleFraction is NaN and its `AF == 'NaN'` rule would no-call a hom-ref genotype.
A sample without data (no read, GT ./.) stays missing. Every fixed ID must come out once.
"""
import argparse
import sys
from collections import defaultdict

import pandas as pd
import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--variants', required=True)
    ap.add_argument('--allsites', required=True)
    ap.add_argument('--fixed', required=True, help='fixed_sites.tsv')
    ap.add_argument('--chrom', help='restrict to one chromosome')
    ap.add_argument('--out', required=True, help='uncompressed VCF')
    a = ap.parse_args()

    fs = pd.read_csv(a.fixed, sep='\t', dtype=str)
    if a.chrom:
        fs = fs[fs.chrom == a.chrom]
    want = defaultdict(list)
    for r in fs.itertuples(index=False):
        want[(r.chrom, int(r.pos))].append((r.variant_id, r.ref, r.alt))

    var = pysam.VariantFile(a.variants)
    alls = pysam.VariantFile(a.allsites)
    samples = list(var.header.samples)
    if list(alls.header.samples) != samples:
        sys.exit('ABORT: the two joint calls hold different samples')

    by_id, by_pos = {}, defaultdict(list)
    for rec in var.fetch():
        if (rec.chrom, rec.pos) not in want:
            continue
        if len(rec.alts or ()) != 1:
            sys.exit(f'ABORT: {rec.id} is not split')
        by_id[rec.id] = rec
        by_pos[(rec.chrom, rec.pos)].append(rec)
    # non-variant templates: an all-sites record in which no sample carries a non-reference allele.
    # Besides ALT '.', GenotypeGVCFs 4.1.4.0 writes such sites with a candidate allele and <NON_REF>
    # (e.g. A -> T,<NON_REF>) when every genotype is 0/0; the variant call then has no record there.
    ref_only = {}
    for rec in alls.fetch():
        key = (rec.chrom, rec.pos)
        if key not in want or key in ref_only:
            continue
        nonref = any(x not in (0, None) for s in samples for x in (rec.samples[s].get('GT') or ()))
        if not nonref:
            ref_only[key] = rec

    hdr = var.header.copy()
    for k, line in (('AD', '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">'),
                    ('GQ', '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'),
                    ('DP', '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">')):
        if k not in hdr.formats:
            hdr.add_line(line)
    hdr.add_line('##platform_depth_expand_fixed=one record per fixed site; 0/0 records built from another '
                 'allele\'s record (GQ, DP, ref AD kept) or from the non-variant record (GQ<-RGQ, AD<-(DP,0))')
    out = pysam.VariantFile(a.out, 'w', header=hdr)
    n = {'kept': 0, 'other_allele': 0, 'non_variant': 0, 'no_record': 0}

    def build(chrom, pos, vid, ref, alt, src, how):
        rec = out.new_record(contig=chrom, start=pos - 1, alleles=(ref, alt), id=vid, qual=None, filter=None)
        for s in samples:
            g = src.samples[s]
            gt = g.get('GT')
            dp = g.get('DP')
            if gt is None or any(x is None for x in gt):
                rec.samples[s]['GT'] = (None, None)
                continue
            rec.samples[s]['GT'] = (0, 0)
            if how == 'other_allele':
                ad = g.get('AD')
                rec.samples[s]['AD'] = (ad[0] if ad and ad[0] is not None else (dp or 0), 0)
                gq = g.get('GQ')
            else:
                ad = g.get('AD') if 'AD' in g else None
                rec.samples[s]['AD'] = (ad[0] if ad and ad[0] is not None else (dp or 0), 0)
                gq = g.get('RGQ') if 'RGQ' in g and g.get('RGQ') is not None else g.get('GQ')
            if dp is not None:
                rec.samples[s]['DP'] = dp
            if gq is not None:
                rec.samples[s]['GQ'] = gq
        return rec

    order = {c: i for i, c in enumerate(hdr.contigs)}
    for (chrom, pos) in sorted(want, key=lambda k: (order[k[0]], k[1])):
        for vid, ref, alt in want[(chrom, pos)]:
            if vid in by_id:
                rec = by_id[vid]
                rec.translate(hdr)
                out.write(rec)
                n['kept'] += 1
            elif by_pos.get((chrom, pos)):
                out.write(build(chrom, pos, vid, ref, alt, by_pos[(chrom, pos)][0], 'other_allele'))
                n['other_allele'] += 1
            elif (chrom, pos) in ref_only:
                out.write(build(chrom, pos, vid, ref, alt, ref_only[(chrom, pos)], 'non_variant'))
                n['non_variant'] += 1
            else:
                n['no_record'] += 1
    out.close()
    print(f'[expand_fixed] {a.chrom or "all"}: {n}')
    if n['no_record']:
        sys.exit(f'ABORT: {n["no_record"]} fixed sites have no record in either joint call')


if __name__ == '__main__':
    main()
