"""summaries of the CRAM scan (all.out.tsv) by cohort and group: ADARB2 carrier-signature reads, LINC01146 base counts"""
import sys
import numpy as np
import pandas as pd
W = sys.argv[1]
T = pd.read_csv('/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/platform_depth/results/00.samples/sample_table.tsv', sep='\t')
o = pd.read_csv(f'{W}/all.out.tsv', sep='\t')
print('scanned', len(o), 'errors', o.error.notna().sum())
ga = pd.read_csv(f'{W}/adarb2.samples.tsv', sep='\t')[['sample_id', 'call', 'ad', 'dp']].rename(columns={'call': 'ad_call', 'ad': 'ad_ad', 'dp': 'ad_dp'})
gl = pd.read_csv(f'{W}/linc.samples.tsv', sep='\t')[['sample_id', 'call', 'ad', 'dp']].rename(columns={'call': 'li_call', 'ad': 'li_ad', 'dp': 'li_dp'})
d = T.merge(o, on='sample_id').merge(ga, on='sample_id').merge(gl, on='sample_id')
d['grp'] = np.where(d.status == 'control', 'controls HiSeq X', 'cases ' + d.platform)
d.to_csv(f'{W}/merged.tsv', sep='\t', index=False)
pd.set_option('display.width', 220)
print(d.groupby('grp').agg(n=('sample_id', 'size'), rlen=('ad_rlen_sum', 'sum'), reads=('ad_reads', 'sum')).assign(rlen=lambda x: x.rlen / x.reads))
# ADARB2: signature reads by called genotype and group
for g, x in d.groupby('grp'):
    print(g, {c: sorted(y.ad_alt.tolist())[-12:] if c != '0/0' else (int((y.ad_alt >= 3).sum()), len(y)) for c, y in x.groupby('ad_call')})
