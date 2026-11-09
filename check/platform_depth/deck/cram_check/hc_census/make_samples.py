#!/usr/bin/env python3
"""The genomes of the HaplotypeCaller census at the ADARB2 lead (user 2026-11-07: "我们有没有什么办法取样取的更加系统
有说服力"): no selection. Every full-mainland genome that carries the deletion in its reads (>= 4 signature reads,
../cram_scan.tsv) or in the cohort's calls (0/1 or 1/1 before QC), and, as negative controls, 20 genomes per platform
drawn at random (fixed seed) from the clear non-carriers (no signature read, called 0/0). Full mainland holds narrow and
intermediate mainland, so the census covers all three cohorts.

Writes samples.tsv (INTERNAL: sample IDs and CRAM paths; git-ignored, never copied into a tracked file) and prints the
counts by group. Usage: python3 make_samples.py"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
MIN_READS = 4                       # as collect_stats.py CRAM_MIN and the pilot (../hc_debug)
N_NONCARRIER = 20                   # random non-carriers per platform
SEED = 20261107
PLATFORMS = ['HiSeqX 15x', 'DNBSeq-T7 30x', 'NovaSeq 30x', 'DNBSeq-G400RS 30x', 'DNBseq-G400RS 15x']

c = pd.read_csv(os.path.join(HERE, '..', 'cram_scan.tsv'), sep='\t', dtype={'sample_id': str})
st = pd.read_csv(os.path.join(HERE, '..', '..', '..', 'results', '00.samples', 'sample_table.tsv'), sep='\t',
                 dtype={'sample_id': str})
c = c[c.in_full_mainland == 1].copy()
assert c.error.isna().all(), 'a scanned genome has an error'
assert set(c.platform) == set(PLATFORMS), f'platforms changed: {sorted(set(c.platform))}'
c['read'] = c.ad_alt >= MIN_READS
c['called'] = c.ad_call.isin(['0/1', '1/1'])
role = pd.Series('', index=c.index)
role[c.read & c.called] = 'read carrier, called'
role[c.read & ~c.called] = 'read carrier, missed'
role[~c.read & c.called] = 'called, no signature'
c['role'] = role
carriers = c[c.role != '']
pool = c[(c.ad_alt == 0) & (c.ad_call == '0/0')]
neg = pd.concat([pool[pool.platform == p].sort_values('sample_id').sample(n=N_NONCARRIER, random_state=SEED)
                 for p in PLATFORMS])
neg = neg.assign(role='non-carrier, random')
d = pd.concat([carriers, neg])
d = d.merge(st[['sample_id', 'cram', 'crai']], on='sample_id', how='left', validate='one_to_one')
assert d.cram.notna().all() and d.crai.notna().all(), 'a genome without its CRAM'
assert d.sample_id.is_unique
d = d[['sample_id', 'role', 'platform', 'status', 'cram', 'crai', 'ad_alt', 'ad_call', 'ad_reads']]
d.to_csv(os.path.join(HERE, 'samples.tsv'), sep='\t', index=False)
print(d.groupby(['platform', 'role']).size().to_string())
print('total', len(d))
