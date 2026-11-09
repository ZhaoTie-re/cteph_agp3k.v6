#!/usr/bin/env python3
"""Stage 01 — read length and depth of a chunk of genomes, measured from their CRAMs.

Depth is what HaplotypeCaller can use: `samtools coverage -q 20 -Q 10` (mapping quality >= 20,
base quality >= 10; note samtools coverage reads -q as MQ and -Q as BQ, the reverse of
samtools depth), unmapped / secondary / QC-fail / duplicate reads excluded. Each region is
queried through the index (`samtools depth -b` streams the whole file instead).

Per genome:
  probe windows  the July audit's ten 1 Mb windows, with duplicates excluded (the effective
                 depth the target is set on) and kept (comparable to the FASTQ-based
                 Observed_Depth);
  targets        every calling interval of stage 02, duplicates excluded;
  read length    of 1,000 primary records (-F 0x900: secondary and supplementary records carry
                 hard-clipped sequence).
"""
import argparse
import statistics
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402

FF_NODUP = 'UNMAP,SECONDARY,QCFAIL,DUP'
FF_DUP = 'UNMAP,SECONDARY,QCFAIL'


def coverage(samtools, cram, crai, region, ff):
    out = subprocess.run([samtools, 'coverage', '-H', '-q', str(K.MIN_MQ), '-Q', str(K.MIN_BQ), '--ff', ff,
                          '--reference', str(K.FASTA), '-r', region, f'{cram}##idx##{crai}'],
                         check=True, capture_output=True, text=True).stdout.strip().split('\n')
    if len(out) != 1:
        raise RuntimeError(f'{cram} {region}: {len(out)} coverage lines')
    f = out[0].split('\t')
    # rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq
    return {'numreads': int(f[3]), 'covbases': int(f[4]), 'meandepth': float(f[6]),
            'meanbaseq': float(f[7]), 'meanmapq': float(f[8])}


def read_length(samtools, cram, crai):
    ch, s, _ = K.PROBE_WINDOWS[0]
    p = subprocess.run(f"{samtools} view -F 0x900 --reference {K.FASTA} '{cram}##idx##{crai}' "
                       f"{ch}:{s}-{s + 200000} | head -n 1000 | cut -f10", shell=True, check=True,
                       capture_output=True, text=True, executable='/bin/bash')
    lens = [len(x) for x in p.stdout.split('\n') if x and x != '*']
    if len(lens) < 100:
        raise RuntimeError(f'{cram}: only {len(lens)} primary reads for the read length')
    return statistics.median(lens), max(lens)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--chunk', required=True, help='TSV: sample_id, cram, crai')
    ap.add_argument('--targets', required=True, help='hc.intervals (chr:start-end per line)')
    ap.add_argument('--samtools', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    chunk = K.read_tsv(a.chunk)
    targets = [ln.strip() for ln in Path(a.targets).read_text().splitlines() if ln.strip()]
    probes = [f'{c}:{s}-{e}' for c, s, e in K.PROBE_WINDOWS]
    with open(a.out, 'w') as fh:
        fh.write('sample_id\tset\tregion\tduplicates\tnumreads\tcovbases\tmeandepth\tmeanbaseq\tmeanmapq\t'
                 'readlen_median\treadlen_max\n')
        for r in chunk.itertuples(index=False):
            med, mx = read_length(a.samtools, r.cram, r.crai)
            jobs = [('probe', reg, 'excluded', FF_NODUP) for reg in probes] + \
                   [('probe', reg, 'kept', FF_DUP) for reg in probes] + \
                   [('target', reg, 'excluded', FF_NODUP) for reg in targets]
            for st, reg, dup, ff in jobs:
                c = coverage(a.samtools, r.cram, r.crai, reg, ff)
                fh.write(f"{r.sample_id}\t{st}\t{reg}\t{dup}\t{c['numreads']}\t{c['covbases']}\t"
                         f"{c['meandepth']:.4f}\t{c['meanbaseq']:.2f}\t{c['meanmapq']:.2f}\t{med}\t{mx}\n")
            fh.flush()
            print(f'[depth] {r.sample_id} done', flush=True)


if __name__ == '__main__':
    main()
