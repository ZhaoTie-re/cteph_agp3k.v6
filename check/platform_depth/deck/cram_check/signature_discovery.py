#!/usr/bin/env python3
"""ADARB2: find the read signature of the carrier haplotype (user 2026-11-06, CRAM check).

The lead chr10:1747167 CGTG>C lies in a minisatellite of about 390 bp, 31-bp units; the aligner places carrier reads without
a gap (one unit shifted, a few mismatches), so the CIGAR does not show the call. Instead, k-mers are ranked by how well
they separate HiSeq X controls called 1/1 and 0/1 from those called 0/0 (reads of the minisatellite, any MAPQ; both
strands). Input: a pickle from dump_reads.py and a table of the dumped samples (sample_id, status, platform, call).
Output (stdout): the best k-mers for K = 15, 19, 21 with the share of samples holding them per group.
Used for the deck: the 16-mer GATCGCATGATGTGTA (ALT) and its reference form GATCGCGTGATGTGTA (three units).
"""
import collections
import pickle
import sys

import pandas as pd

R = pickle.load(open(sys.argv[1], 'rb'))
meta = pd.read_csv(sys.argv[2], sep='\t').set_index('sample_id')
comp = str.maketrans('ACGTN', 'TGCAN')


def canon(k):
    rc = k.translate(comp)[::-1]
    return min(k, rc)


grp = meta.apply(lambda r: f"{r.status}|{'T7' if 'T7' in r.platform else ('HiSeqX' if r.status == 'control' else 'otherCase')}"
                           f"|{r.call}", axis=1)
G = {g: list(grp[grp == g].index) for g in grp.unique()}
for K in (15, 19, 21):
    KS = {s: {canon(seq[i:i + K]) for *_, seq in R[s] if seq for i in range(len(seq) - K + 1)} for s in R}
    allk = collections.Counter(k for s in R for k in KS[s])

    def frac(k, g):
        return sum(k in KS[s] for s in G.get(g, [])) / max(1, len(G.get(g, [])))
    rows = [(k, frac(k, 'control|HiSeqX|1/1'), frac(k, 'control|HiSeqX|0/1'), frac(k, 'control|HiSeqX|0/0'),
             frac(k, 'case|otherCase|0/1'), frac(k, 'case|T7|0/0')) for k, n in allk.items() if n >= 5]
    d = pd.DataFrame(rows, columns=['kmer', 'hom', 'het', 'ref', 'other_case_carriers', 't7'])
    d['score'] = (d.hom + d.het) / 2 - d.ref
    print('K =', K)
    print(d.sort_values('score', ascending=False).head(12).round(2).to_string(index=False))
