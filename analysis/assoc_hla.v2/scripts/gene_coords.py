#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : GRCh38 span of every typed HLA gene, from the NCBI RefSeq gene table
#           ANNOVAR ships (hg38_refGene.txt, UCSC genePred with a leading bin).
#
#           The span is the union of the gene's transcripts on the primary chr6
#           (alternate-haplotype contigs such as chr6_GL000250v2_alt are
#           ignored). It is used for two things only: to order the genes along
#           the chromosome in tables and figures, and to give each marker a
#           unique chromosome-6 position in the PLINK sets (gene start + rank of
#           the marker within the gene). No association result depends on it.
#
#           Every gene in --genes must be found exactly once on chr6; a missing
#           gene stops the run rather than placing its markers nowhere.
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import csv
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--refgene', required=True, help='hg38_refGene.txt (genePred with bin column)')
    p.add_argument('--genes', required=True, help='comma-separated IPD-IMGT/HLA locus names')
    p.add_argument('--out', required=True)
    return p.parse_args()


def main():
    a = parse_args()
    genes = [g.strip() for g in a.genes.split(',') if g.strip()]
    want = {f'HLA-{g}': g for g in genes}
    span = {}
    with open(a.refgene) as fh:
        for row in csv.reader(fh, delimiter='\t'):
            # bin, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd,
            # exonCount, exonStarts, exonEnds, score, name2, ...
            if row[2] != 'chr6' or row[12] not in want:
                continue
            g = want[row[12]]
            s, e = int(row[4]), int(row[5])
            if g in span:
                s0, e0, st0, n0 = span[g]
                if st0 != row[3]:
                    sys.exit(f'ERROR: {row[12]} has transcripts on both strands of chr6')
                span[g] = (min(s0, s), max(e0, e), st0, n0 + 1)
            else:
                span[g] = (s, e, row[3], 1)
    missing = [g for g in genes if g not in span]
    if missing:
        sys.exit(f'ERROR: no chr6 RefSeq transcript for HLA-{", HLA-".join(missing)} in {a.refgene}')

    order = sorted(genes, key=lambda g: span[g][0])
    with open(a.out, 'w') as out:
        out.write('gene\tsymbol\tchrom\tgene_start\tgene_end\tstrand\ttranscripts\torder\n')
        for i, g in enumerate(order, 1):
            s, e, st, n = span[g]
            # genePred txStart is 0-based; written here 1-based, closed.
            out.write(f'{g}\tHLA-{g}\t6\t{s + 1}\t{e}\t{st}\t{n}\t{i}\n')
    print(f'[gene_coords] {len(order)} genes: ' + ' '.join(order))


if __name__ == '__main__':
    main()
