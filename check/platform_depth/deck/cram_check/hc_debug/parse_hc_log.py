#!/usr/bin/env python3
"""Summarize the HaplotypeCaller debug log (--debug-assembly) of each genome at the ADARB2 lead: every assembly region
overlapping the lead window, its reads, the k-mer sizes tried and why each was dropped, the candidate haplotypes (and
whether any carries the one-unit deletion: a CIGAR with D), the events found and the trimmed span.
Usage: python3 parse_hc_log.py > hc_summary.tsv"""
import re
import sys
import pandas as pd

LEAD = (1747150, 1747210)         # the lead (1747167) and the two units holding the G>A of the signature (1747168, 1747199)
S = pd.read_csv('samples.tsv', sep='\t', dtype={'sample_id': str})
rows = []
for _, s in S.iterrows():
    blocks, cur = [], None
    for line in open(f'out/{s.sample_id}/{s.sample_id}.hc.log'):
        m = re.search(r'Assembling (chr10):(\d+)-(\d+) with (\d+) reads', line)
        if m:
            cur = {'start': int(m.group(2)), 'end': int(m.group(3)), 'reads': int(m.group(4)), 'used': [], 'dropped': [],
                   'n_hap': None, 'del_haps': 0, 'cigars': set(), 'events': [], 'trim': ''}
            blocks.append(cur)
            continue
        if cur is None:
            continue
        if (m := re.search(r'Using kmer size of (\d+)', line)):
            cur['used'].append(int(m.group(1)))
        elif (m := re.search(r'Not using kmer size of (\d+) in read threading assembler because (.*)', line)):
            cur['dropped'].append(f'{m.group(1)}:{"non-unique" if "non-unique" in m.group(2) else m.group(2).strip()}')
        elif (m := re.search(r'Found (\d+) candidate haplotypes', line)):
            cur['n_hap'] = int(m.group(1))
        elif (m := re.search(r'> Cigar = (\S+) :', line)):
            cur['cigars'].add(m.group(1))
            cur['del_haps'] += 'D' in m.group(1)
        elif (m := re.search(r'Genotyping event at (\d+) with alleles = (\[.*\])', line)):
            cur['events'].append(f'{m.group(1)}{m.group(2).replace(" ", "")}')
        elif (m := re.search(r'Trimmed region to (chr10:\S+)', line)):
            cur['trim'] = m.group(1)
    for b in blocks:
        if b['start'] <= LEAD[1] and b['end'] >= LEAD[0]:
            rows.append({'sample_id': s.sample_id, 'role': s.role, 'region': f"{b['start']}-{b['end']}", 'reads': b['reads'],
                         'kmers_used': ','.join(map(str, b['used'])), 'kmers_dropped': ' '.join(b['dropped']),
                         'haplotypes': b['n_hap'], 'haps_with_deletion': b['del_haps'],
                         'cigars': ' '.join(sorted(b['cigars'])), 'events': ' '.join(b['events']), 'trimmed_to': b['trim']})
pd.DataFrame(rows).to_csv(sys.stdout, sep='\t', index=False)
