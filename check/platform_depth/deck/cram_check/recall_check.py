#!/usr/bin/env python3
"""Where the ADARB2 deletion is lost (deck, user 2026-11-06: "我们有没有办法确认是calling的问题？", then "26～29 重构").
Read only. Two tables for collect_stats.py (stats.json cram.aligner, cram.recall):

1. aligner.tsv: the alignment program records (@PG) of every full-mainland CRAM. The bwa command is normalized (program
   and reference by base name; the thread count and the read group, which name the sample, dropped), so one command
   for every genome means one alignment for every platform.
2. recall_adarb2.tsv: the genotype at the ADARB2 lead (chr10:1747167:CGTG:C) after Q1's re-calling of the deep cases'
   reads at full depth (the baseline arm, f = 1: every read, GATK HaplotypeCaller 4.1.4.0 with the cohort's settings,
   joint genotyping, the cohort's QC; results/03.calls/downsampling/arms.qc.vcf.gz), with the CRAM scan's
   signature-read count and the cohort's original call (cram_scan.tsv). A reference-only record becomes 0/0 in the
   arms (scripts/expand_fixed.py), so 0/0 means the caller wrote no deletion allele for that genome.

Usage: python3 -B recall_check.py [n_proc]
"""
import os
import subprocess
import sys
from multiprocessing import Pool

import pandas as pd
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SAMPLES = os.path.join(ROOT, 'results', '00.samples', 'sample_table.tsv')
ARMS = os.path.join(ROOT, 'results', '03.calls', 'downsampling', 'arms.qc.vcf.gz')
LEAD = ('chr10', 1747167, 'chr10:1747167:CGTG:C')
BCFTOOLS = '/home/b/b37974/bcftools/bcftools'


def norm_bwa(cl):
    out, skip = [], False
    for t in cl.split():
        if skip:
            skip = False
            continue
        if t.startswith('-t'):                    # thread count
            if t == '-t':
                skip = True
            continue
        if t.startswith('-R'):                    # read group: names the sample
            if t == '-R':
                skip = True
            continue
        out.append(os.path.basename(t) if '/' in t else t)
    return ' '.join(out)


def header(args):
    sid, cram = args
    try:
        h = pysam.AlignmentFile(cram, 'rc').header.to_dict()
    except Exception as e:  # noqa: BLE001
        return {'sample_id': sid, 'error': str(e)}
    pg = h.get('PG', [])
    bwa = [p for p in pg if p.get('ID', '').startswith('bwa') or p.get('PN') == 'bwa']
    sb = [p for p in pg if p.get('ID', '').upper().startswith('SAMBLASTER')]
    return {'sample_id': sid, 'error': '', 'n_bwa': len(bwa),
            'bwa_vn': bwa[0].get('VN', '') if bwa else '', 'bwa_cmd': norm_bwa(bwa[0].get('CL', '')) if bwa else '',
            'samblaster_vn': sb[0].get('VN', '') if sb else '',
            'pl': ';'.join(sorted({r.get('PL', '') for r in h.get('RG', [])}))}


def main():
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    st = pd.read_csv(SAMPLES, sep='\t', dtype={'sample_id': str})
    st = st[st.in_full_mainland.astype(str) == '1']
    with Pool(nproc) as pool:
        rows = pool.map(header, list(zip(st.sample_id, st.cram)), chunksize=8)
    a = st[['sample_id', 'status', 'platform']].merge(pd.DataFrame(rows), on='sample_id')
    a.to_csv(os.path.join(HERE, 'aligner.tsv'), sep='\t', index=False)
    print(f'aligner: {len(a)} genomes, {int((a.error != "").sum())} errors')
    print(a.groupby(['platform', 'bwa_vn', 'bwa_cmd', 'samblaster_vn']).size().to_string())
    # the re-called genotype at the lead, baseline arm
    q = subprocess.run([BCFTOOLS, 'view', '-t', f'{LEAD[0]}:{LEAD[1]}', '-i', f'ID="{LEAD[2]}"', ARMS],
                       check=True, capture_output=True)
    g = subprocess.run([BCFTOOLS, 'query', '-f', '[%SAMPLE\t%GT\n]', '-'], input=q.stdout, check=True,
                       capture_output=True).stdout.decode().splitlines()
    gt = pd.DataFrame([x.split('\t') for x in g], columns=['arm_sample', 'gt_recall'])
    gt = gt[gt.arm_sample.str.endswith('__baseline')].assign(sample_id=lambda d: d.arm_sample.str.replace('__baseline', '', regex=False))
    if gt.empty:
        sys.exit('ABORT: no baseline genotype at the ADARB2 lead')
    cs = pd.read_csv(os.path.join(HERE, 'cram_scan.tsv'), sep='\t', dtype={'sample_id': str})
    r = cs[['sample_id', 'status', 'platform', 'ad_alt', 'ad_call']].merge(gt[['sample_id', 'gt_recall']], on='sample_id')
    r.to_csv(os.path.join(HERE, 'recall_adarb2.tsv'), sep='\t', index=False)
    print(f'recall: {len(r)} genomes re-called at full depth')
    print(r.assign(read_carrier=r.ad_alt >= 4).groupby(['platform', 'read_carrier', 'ad_call', 'gt_recall']).size().to_string())


if __name__ == '__main__':
    main()
