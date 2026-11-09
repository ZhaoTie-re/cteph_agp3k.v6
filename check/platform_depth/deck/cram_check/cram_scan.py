#!/usr/bin/env python3
"""Read-level evidence at the two flagged leads, straight from each sample's CRAM (read only).

ADARB2 chr10:1747167 CGTG>C (GTG deleted at 1-based 1747168-1747170) lies in a minisatellite of about 390 bp, 31-bp units;
the same 3-bp deletion can be placed in another unit (31 bp apart), so every 3-bp deletion inside the repeat is
counted, with its position. LINC01146 chr14:88019710 G>A: base counts by strand and quality.
Reads: primary, mapped, not duplicate, not QC-fail; MAPQ and base quality kept so filters can be applied later.
Output: one TSV row per sample and locus.
"""
import sys
import time
from multiprocessing import Pool

import pysam

FA = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa'
AD_CHR, AD_POS0 = 'chr10', 1747167            # 0-based first deleted base (1-based 1747168)
VNTR = (1746850, 1747320)                     # 0-based half-open window of the minisatellite and flanks
ALT, REF = 'GATCGCATGATGTGTA', 'GATCGCGTGATGTGTA'
_C = str.maketrans('ACGTN', 'TGCAN')
ALT_RC, REF_RC = ALT.translate(_C)[::-1], REF.translate(_C)[::-1]
LI_CHR, LI_POS0 = 'chr14', 88019709           # 0-based position of 1-based 88019710


def scan(args):
    sid, cram, crai = args
    t0 = time.time()
    out = {'sample_id': sid}
    try:
        f = pysam.AlignmentFile(cram, 'rc', reference_filename=FA, index_filename=crai)
    except Exception as e:  # noqa: BLE001
        return {**out, 'error': f'open: {e}'}
    # ADARB2: alignment-free. The carrier haplotype is read by its signature (the unit "GATCGCATGATGTGTA", G>A in a
    # perfect 31-bp unit; its reference form "GATCGCGTGATGTGTA" sits in three units), counted in every read of the
    # minisatellite whatever its alignment (any MAPQ; unmapped mates placed here included)
    n = dict(ad_reads=0, ad_rlen_sum=0, ad_alt=0, ad_alt_q20=0, ad_alt_clip=0, ad_ref=0, ad_ref_q20=0, ad_mapq_lt20=0)
    for r in f.fetch(AD_CHR, VNTR[0], VNTR[1]):
        if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail:
            continue
        seq = r.query_sequence or ''
        n['ad_reads'] += 1
        n['ad_rlen_sum'] += len(seq)
        q20 = (not r.is_unmapped) and r.mapping_quality >= 20
        n['ad_mapq_lt20'] += not q20
        if ALT in seq or ALT_RC in seq:
            n['ad_alt'] += 1
            n['ad_alt_q20'] += q20
            n['ad_alt_clip'] += 'S' in (r.cigarstring or '')
        if REF in seq or REF_RC in seq:
            n['ad_ref'] += 1
            n['ad_ref_q20'] += q20
    out.update(n)
    # LINC01146
    c = dict(li_reads=0, li_ref=0, li_alt=0, li_other=0, li_ref_q=0, li_alt_q=0, li_alt_fwd_q=0, li_ref_fwd_q=0,
             li_mapq_lt20=0)
    for r in f.fetch(LI_CHR, LI_POS0, LI_POS0 + 1):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail:
            continue
        qp = None
        for qpos, rpos in r.get_aligned_pairs(matches_only=True):
            if rpos == LI_POS0:
                qp = qpos
                break
        if qp is None:
            continue
        c['li_reads'] += 1
        base, bq = r.query_sequence[qp], r.query_qualities[qp]
        good = r.mapping_quality >= 20 and bq >= 20
        c['li_mapq_lt20'] += r.mapping_quality < 20
        if base == 'G':
            c['li_ref'] += 1
            c['li_ref_q'] += good
            c['li_ref_fwd_q'] += good and not r.is_reverse
        elif base == 'A':
            c['li_alt'] += 1
            c['li_alt_q'] += good
            c['li_alt_fwd_q'] += good and not r.is_reverse
        else:
            c['li_other'] += 1
    out.update(c)
    out['seconds'] = round(time.time() - t0, 2)
    f.close()
    return out


def main():
    lst, outp, nproc = sys.argv[1], sys.argv[2], int(sys.argv[3])
    jobs = [tuple(l.rstrip('\n').split('\t')[:3]) for l in open(lst) if l.strip()]
    cols = None
    with Pool(nproc) as pool, open(outp, 'w') as fo:
        for i, row in enumerate(pool.imap_unordered(scan, jobs, chunksize=1)):
            if cols is None:
                cols = ['sample_id', 'error', 'ad_reads', 'ad_rlen_sum', 'ad_alt', 'ad_alt_q20', 'ad_alt_clip', 'ad_ref',
                        'ad_ref_q20', 'ad_mapq_lt20', 'li_reads', 'li_ref', 'li_alt', 'li_other', 'li_ref_q', 'li_alt_q',
                        'li_ref_fwd_q', 'li_alt_fwd_q', 'li_mapq_lt20', 'seconds']
                fo.write('\t'.join(cols) + '\n')
            fo.write('\t'.join(str(row.get(k, '')) for k in cols) + '\n')
            fo.flush()


if __name__ == '__main__':
    main()
