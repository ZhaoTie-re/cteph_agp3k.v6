#!/usr/bin/env python3
"""LINC01146 and the G400RS 15x cases (user 2026-11-06: "LINC01146 应该是G400RS 15X的问题，你看看"). Read only.

1. linc_fragments.tsv: at chr14:88019710, per case genome (all case platforms), every primary mapped read covering the
   site: reads, duplicate-flagged reads, independent fragments (start, mate start, strand) and G / A among them; a
   homozygous call is allelic dropout only if few independent fragments carry it.
2. het_by_sample.tsv: genome-wide heterozygosity and genotype missingness of every full-mainland genome on the
   pipeline's common null SNVs (results/02.targets/null_snps.tsv, SNVs with MAF bin >= 0.1: 6,000), from the cohort's
   QC'd genotypes (plink2 --het --missing; frequencies from all 3,101 genomes).
Usage: python3 linc_g400rs15x.py <out_dir>
"""
import os
import subprocess
import sys

import pandas as pd
import pysam

OUT = sys.argv[1]
COMP = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/platform_depth'
FA = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa'
BASE = ('/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/12_model_inputs/full_mainland/'
        'genotype/base/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm')
PLINK2 = '/home/b/b37974/plink2_alpha6/plink2'
T = pd.read_csv(f'{COMP}/results/00.samples/sample_table.tsv', sep='\t', dtype={'sample_id': str})
T = T[T.in_full_mainland == 1]

rows = []
for r in T[T.status == 'case'].itertuples():
    f = pysam.AlignmentFile(r.cram, 'rc', reference_filename=FA, index_filename=r.crai)
    reads = dups = 0
    frags, base = set(), {}
    for a in f.fetch('chr14', 88019709, 88019710):
        if a.is_secondary or a.is_supplementary or a.is_unmapped:
            continue
        q = [qp for qp, rp in a.get_aligned_pairs(matches_only=True) if rp == 88019709]
        if not q:
            continue
        reads += 1
        if a.is_duplicate:
            dups += 1
            continue
        frags.add((a.reference_start, a.next_reference_start, a.is_reverse))
        b = a.query_sequence[q[0]]
        base[b] = base.get(b, 0) + 1
    f.close()
    rows.append((r.sample_id, r.platform, reads, dups, len(frags), base.get('G', 0), base.get('A', 0)))
pd.DataFrame(rows, columns=['sample_id', 'platform', 'reads', 'dup_flagged', 'fragments', 'G', 'A']).to_csv(
    os.path.join(OUT, 'linc_fragments.tsv'), sep='\t', index=False)

nul = pd.read_csv(f'{COMP}/results/02.targets/null_snps.tsv', sep='\t')
sel = nul[(~nul.is_indel) & nul.maf_bin.isin(['0.1-0.2', '0.2-0.5'])].variant_id
tmp = os.path.join(OUT, 'null_snv_ge10.txt')
sel.to_csv(tmp, index=False, header=False)
pre = os.path.join(OUT, 'het_null')
subprocess.run([PLINK2, '--bfile', BASE, '--extract', tmp, '--het', '--missing', '--memory', '6000', '--threads', '2',
                '--out', pre], check=True, stdout=subprocess.DEVNULL)
h = pd.read_csv(pre + '.het', sep='\t', dtype={'IID': str, '#FID': str}).rename(columns={'IID': 'sample_id'})
m = pd.read_csv(pre + '.smiss', sep='\t', dtype={'IID': str, '#FID': str}).rename(columns={'IID': 'sample_id'})
d = h.merge(m[['sample_id', 'F_MISS']], on='sample_id').merge(
    T[['sample_id', 'status', 'platform', 'observed_depth']], on='sample_id', how='left')
assert d.platform.notna().all() and len(d) == len(T)
d[['sample_id', 'status', 'platform', 'observed_depth', 'O(HOM)', 'E(HOM)', 'OBS_CT', 'F', 'F_MISS']].to_csv(
    os.path.join(OUT, 'het_by_sample.tsv'), sep='\t', index=False)
for f in ('.het', '.smiss', '.vmiss', '.log'):
    if os.path.exists(pre + f):
        os.remove(pre + f)
os.remove(tmp)
print('wrote linc_fragments.tsv and het_by_sample.tsv in', OUT)
