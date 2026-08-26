#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The list of reference sequences an HLA read can be sitting on, built
#           once for the whole run from the .fai index.
#
#           Three groups plus one wildcard, and every one of them is load-bearing:
#             - the classical MHC window on chr6
#             - the chr6 ALT scaffolds, which carry divergent haplotypes
#             - the 525 HLA* decoy contigs of hs38DH; reads from an allele far
#               from the primary reference map HERE and nowhere else
#             - '*', the unmapped reads, which is where a read whose allele is in
#               none of the above ends up
#           Dropping any one of them loses reads from exactly the samples whose
#           alleles are most unusual, which is the opposite of what typing needs.
#
#           Read from hs38DH.fa.fai (151 KB), not by parsing the 3 GB FASTA.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fai', required=True, help='<reference>.fa.fai')
    p.add_argument('--mhc-chrom', default='chr6')
    p.add_argument('--mhc-start', type=int, default=29540000)
    p.add_argument('--mhc-end', type=int, default=33420000)
    p.add_argument('--out', default='contig_list.txt')
    return p.parse_args()


def main():
    args = parse_args()
    fai = Path(args.fai)
    if not fai.exists():
        raise SystemExit(f'ABORT: no .fai at {fai}. samtools faidx the reference first.')

    names = [ln.split('\t', 1)[0] for ln in fai.read_text().splitlines() if ln.strip()]
    if not names:
        raise SystemExit(f'ABORT: {fai} is empty')

    mhc = f'{args.mhc_chrom}:{args.mhc_start}-{args.mhc_end}'
    alt = [n for n in names if n.startswith(args.mhc_chrom) and n != args.mhc_chrom]
    hla = [n for n in names if n.startswith('HLA')]

    if not hla:
        raise SystemExit(
            f'ABORT: {fai} lists no HLA* contig.\n'
            f'  This reference is not HLA-aware (hs38DH and friends carry 525 of them).\n'
            f'  Typing against it would silently lose the reads that matter most.')
    if args.mhc_chrom not in names:
        raise SystemExit(f'ABORT: {fai} has no sequence named {args.mhc_chrom!r}')

    regions = [mhc] + alt + hla + ['*']
    Path(args.out).write_text('\n'.join(regions) + '\n')

    print(f'[contig_list] {len(regions)} region(s) -> {args.out}')
    print(f'    MHC window        {mhc}')
    print(f'    {args.mhc_chrom} ALT scaffolds  {len(alt)}')
    print(f'    HLA* contigs      {len(hla)}')
    print(f'    unmapped          *')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in contig_list: {e}', file=sys.stderr)
        sys.exit(1)
