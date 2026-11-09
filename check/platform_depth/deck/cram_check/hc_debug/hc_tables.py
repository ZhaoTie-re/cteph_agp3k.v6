#!/usr/bin/env python3
"""Tables behind the HaplotypeCaller diagnostics at the ADARB2 lead (README.md). Read only; run after run_hc_debug.sh
and parse_hc_log.py.

1. hc_genotypes.tsv (internal: sample IDs): each diagnosed genome of samples.tsv, whether HaplotypeCaller assembled a
   haplotype carrying the one-unit deletion (hc_summary.tsv: a candidate haplotype whose CIGAR holds a D), the k-mer
   size it assembled with, and the genotype it wrote at the lead (1747167), at the 28-bp deletion (1747171) and at the
   two G>A positions of the signature (1747168, 1747199) in its debug GVCF.
2. cohort_records.tsv: the cohort's records in the lead window before the site QC (wgs.auto.par 01_prepare_vcf, VQSR
   PASS) and whether each survives it (02_filter_vqc: VQSLOD > 10 & MQ > 58.75), with VQSLOD, MQ, AC, AN.
3. cohort_ga.tsv (counts only): the full-mainland genomes of ../cram_scan.tsv by group, their read carriers
   (>= 4 signature reads), called and missed, and how many of each carry chr10:1747168 G>A in the cohort's calls
   before the site QC.

Usage: python3 hc_tables.py      (writes the three tables next to this script)"""
import os
import pandas as pd
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
V6 = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..', '..'))          # cteph_agp3k.v6
V01 = os.path.join(V6, 'wgs.auto.par/results/01_prepare_vcf/chr10.selected.pass.norm_split.setid.vcf.gz')
V02 = os.path.join(V6, 'wgs.auto.par/results/02_filter_vqc/chr10.vqc.vcf.gz')
CHROM, WIN = 'chr10', (1747160, 1747205)
SITES = {1747167: 'lead', 1747171: 'del28', 1747168: 'ga_1747168', 1747199: 'ga_1747199'}
MIN_READS = 4                                                                     # as collect_stats.py CRAM_MIN


def gt_str(s):
    a = s.get('GT')
    return '/'.join('.' if x is None else str(x) for x in a) if a else './.'


# 1. per diagnosed genome
S = pd.read_csv(os.path.join(HERE, 'samples.tsv'), sep='\t', dtype={'sample_id': str})
H = pd.read_csv(os.path.join(HERE, 'hc_summary.tsv'), sep='\t', dtype={'sample_id': str})
rows = []
for _, s in S.iterrows():
    h = H[H.sample_id == s.sample_id]
    r = {'sample_id': s.sample_id, 'role': s.role, 'platform': s.platform, 'signature_reads': int(s.ad_alt),
         'cohort_call': s.ad_call, 'regions': len(h), 'reads': ';'.join(map(str, h.reads)),
         'kmer_used': ';'.join(map(str, h.kmers_used)), 'deletion_haplotypes': int(h.haps_with_deletion.sum())}
    with pysam.VariantFile(os.path.join(HERE, 'out', s.sample_id, f'{s.sample_id}.g.vcf.gz')) as vf:
        smp = list(vf.header.samples)[0]
        for pos, name in SITES.items():
            hit = [x for x in vf.fetch(CHROM, pos - 1, pos) if x.pos == pos and any(a != '<NON_REF>' for a in x.alts)]
            if hit:
                x = hit[0]
                sx = x.samples[smp]
                r[f'{name}_alleles'] = f"{x.ref}>{','.join(x.alts)}"
                r[f'{name}_gt'] = gt_str(sx)
                r[f'{name}_ad'] = ','.join(map(str, sx.get('AD') or []))
                r[f'{name}_gq'] = sx.get('GQ')
            else:
                r[f'{name}_alleles'], r[f'{name}_gt'], r[f'{name}_ad'], r[f'{name}_gq'] = '', 'ref block', '', ''
    rows.append(r)
G = pd.DataFrame(rows)
G.to_csv(os.path.join(HERE, 'hc_genotypes.tsv'), sep='\t', index=False)

# 2. the cohort's records in the window, before and after the site QC
kept = set()
with pysam.VariantFile(V02) as vf:
    for x in vf.fetch(CHROM, WIN[0] - 1, WIN[1]):
        kept.add((x.pos, x.ref, x.alts[0]))
rec = []
with pysam.VariantFile(V01, drop_samples=True) as vf:
    for x in vf.fetch(CHROM, WIN[0] - 1, WIN[1]):
        i = x.info
        rec.append({'pos': x.pos, 'ref': x.ref, 'alt': x.alts[0], 'filter_01': ';'.join(x.filter.keys()),
                    'vqslod': i.get('VQSLOD'), 'mq': i.get('MQ'),
                    'ac': (i.get('AC') or [None])[0], 'an': i.get('AN'),
                    'kept_after_site_qc': int((x.pos, x.ref, x.alts[0]) in kept)})
pd.DataFrame(rec).to_csv(os.path.join(HERE, 'cohort_records.tsv'), sep='\t', index=False)

# 3. the G>A at 1747168 in the cohort's calls before the site QC, by group (counts only)
CR = pd.read_csv(os.path.join(HERE, '..', 'cram_scan.tsv'), sep='\t', dtype={'sample_id': str})
CR = CR[CR.in_full_mainland == 1].copy()
with pysam.VariantFile(V01) as vf:
    x = [x for x in vf.fetch(CHROM, 1747167, 1747168) if x.pos == 1747168 and x.ref == 'G' and x.alts == ('A',)]
    assert len(x) == 1, 'cohort: not one G>A record at 1747168'
    ga = {sid: gt_str(x[0].samples[sid]) for sid in x[0].samples}
CR['ga'] = CR.sample_id.map(ga)
assert CR.ga.notna().all(), 'cohort: a scanned genome is missing from the cohort calls'
CR['ga_carrier'] = CR.ga.isin(['0/1', '1/1'])
CR['called'] = CR.ad_call.isin(['0/1', '1/1'])
CR['read'] = CR.ad_alt >= MIN_READS
GRP = {'HiSeq X controls': CR.status == 'control', 'DNBSeq-T7 cases': CR.platform == 'DNBSeq-T7 30x',
       'NovaSeq cases': CR.platform == 'NovaSeq 30x', 'G400RS 30x cases': CR.platform == 'DNBSeq-G400RS 30x',
       'G400RS 15x cases': CR.platform == 'DNBseq-G400RS 15x', 'all': CR.status.notna()}
out = []
for g, m in GRP.items():
    y = CR[m]
    miss = y[y.read & ~y.called]
    out.append({'group': g, 'n': len(y), 'read_carriers': int(y.read.sum()), 'called_carriers': int(y.called.sum()),
                'read_carriers_missed': len(miss), 'missed_with_ga': int(miss.ga_carrier.sum()),
                'called_with_ga': int(y[y.called].ga_carrier.sum()),
                'noncarriers_with_ga': int(y[~y.read & ~y.called].ga_carrier.sum()),
                'ga_carriers': int(y.ga_carrier.sum())})
pd.DataFrame(out).to_csv(os.path.join(HERE, 'cohort_ga.tsv'), sep='\t', index=False)
print(G.drop(columns=['reads']).to_string())
print(pd.DataFrame(rec).to_string())
print(pd.DataFrame(out).to_string())
