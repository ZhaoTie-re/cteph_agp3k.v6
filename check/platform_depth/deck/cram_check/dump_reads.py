#!/usr/bin/env python3
"""dump every primary, non-duplicate read (any MAPQ, unmapped mates placed there included) in a region, per sample"""
import pickle
import sys
from multiprocessing import Pool

import pysam

FA = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa'
CHR, A, B = sys.argv[3], int(sys.argv[4]), int(sys.argv[5])


def one(args):
    sid, cram, crai = args
    f = pysam.AlignmentFile(cram, 'rc', reference_filename=FA, index_filename=crai)
    rows = []
    for r in f.fetch(CHR, A, B):
        if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail:
            continue
        rows.append((r.query_name, r.flag, r.reference_start, r.mapping_quality, r.cigarstring, r.query_sequence))
    f.close()
    return sid, rows


if __name__ == '__main__':
    jobs = [tuple(l.rstrip('\n').split('\t')[:3]) for l in open(sys.argv[1]) if l.strip()]
    with Pool(6) as p:
        res = dict(p.map(one, jobs))
    pickle.dump(res, open(sys.argv[2], 'wb'))
    print(len(res), 'samples', sum(len(v) for v in res.values()), 'reads')
