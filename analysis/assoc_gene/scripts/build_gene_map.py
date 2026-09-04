#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The canonical variant -> gene map for ONE cohort x ONE chromosome.
#
#           This is the keystone of the component. Everything downstream reads
#           it: the rvtest --setFile, the per-gene allele and carrier counts,
#           the Bonferroni denominator, and the annotation-distribution figure.
#           Because the set file, the counts and the denominator are all derived
#           from this one table, "the tested sets are the counted sets" becomes
#           a property that can be asserted rather than hoped for.
#
#           HOW THE JOIN WORKS. The callset .bim carries IDs of the exact form
#           chr<N>:<POS>:<REF>:<ALT>, and the JHRP snpEff index is a
#           position-sorted TSV keyed on (#CHROM, POS, REF, ALT) with CHROM
#           already 'chr'-prefixed. So paste(CHROM,':',POS,':',REF,':',ALT)
#           reproduces the bim ID exactly and the join needs NO format
#           conversion. Both inputs are position-sorted, so this is a streaming
#           merge with O(1) memory rather than a 25M-row dictionary.
#
#           The retired VCF path reached the same annotations by exporting
#           PLINK to VCF, running an 8-core bcftools annotate against this same
#           index, and reading the INFO field back. That step was the
#           wall-clock driver of the whole component (203 min). It exists only
#           to move these columns through a VCF, which this component does not
#           need.
#
#           GENE ASSIGNMENT. 29% of variants carry more than one annotation row
#           and those rows may name different genes at different impacts (a
#           variant can be HIGH in gene A and MODERATE in gene B). The default
#           is therefore PER-GENE most-severe: one row per (variant, gene) pair
#           carrying that pair's most severe consequence. A variant genuinely
#           affecting two genes is tested in both, which is what the field does
#           and what rvtest supports -- the same 1-bp range simply appears in
#           two set lines.
# Component: assoc_gene
# ---------------------------------------------------------------------------
import argparse
import gzip
import sys
from collections import Counter, defaultdict

# snpEff's own severity ladder. MODIFIER is last and is never a tested stratum:
# it is what snpEff gives intergenic and deep-intronic variants, which a
# gene-based test has no principled way to assign.
IMPACT_RANK = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}
UNANNOTATED = 'UNANNOTATED'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bim', required=True, help='single-chromosome slice of the callset .bim')
    p.add_argument('--keep-variants', required=True,
                   help='variant IDs surviving --mac (plink2 --write-snplist output)')
    p.add_argument('--snpeff-index', required=True, help='all.VQSR3.chr<N>.vcf_out.tsv.2.gz')
    p.add_argument('--chrom', required=True, help='plain chromosome number, as in bim column 1')
    p.add_argument('--impacts', nargs='+', default=['HIGH', 'MODERATE', 'LOW'],
                   help='impact classes eligible for the map (MODIFIER is never eligible)')
    p.add_argument('--gene-assignment', choices=['per_gene_most_severe', 'variant_most_severe'],
                   default='per_gene_most_severe')
    p.add_argument('--out-map', default=None)
    p.add_argument('--out-anno', default=None)
    p.add_argument('--out-qc', default=None)
    return p.parse_args()


def opener(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def read_bim(path, chrom, keep):
    """Yield (pos, variant_id, ref, alt) in file order; assert the ID convention."""
    prev = -1
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            f = line.split()
            if len(f) < 6:
                raise SystemExit(f'ABORT: {path}:{ln} has {len(f)} fields, expected 6')
            c, vid, pos, a1, a2 = f[0], f[1], int(f[3]), f[4], f[5]
            if c != chrom:
                raise SystemExit(f'ABORT: {path}:{ln} is chromosome {c!r}, expected {chrom!r}')
            # The ID is the join key. Never trust it silently.
            parts = vid.split(':')
            if len(parts) != 4:
                raise SystemExit(f'ABORT: {path}:{ln} id {vid!r} is not chr:pos:ref:alt')
            if parts[0] != f'chr{chrom}' or int(parts[1]) != pos:
                raise SystemExit(f'ABORT: {path}:{ln} id {vid!r} disagrees with '
                                 f'CHROM={c} POS={pos}')
            if {parts[2], parts[3]} != {a1, a2}:
                raise SystemExit(f'ABORT: {path}:{ln} id {vid!r} alleles disagree with '
                                 f'bim A1={a1} A2={a2}')
            if pos < prev:
                raise SystemExit(f'ABORT: {path} is not position-sorted at line {ln} '
                                 f'({pos} < {prev}); the merge join requires sorted input')
            prev = pos
            if keep is not None and vid not in keep:
                continue
            yield pos, vid, parts[2], parts[3]


def index_blocks(path, chrom):
    """Yield (pos, {(ref, alt): [(effect, impact, gene, geneid, biotype), ...]})."""
    want = f'chr{chrom}'
    with opener(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        col = {name: i for i, name in enumerate(head)}
        for need in ('#CHROM', 'POS', 'REF', 'ALT', 'effect', 'impact', 'gene', 'geneid', 'biotype'):
            if need not in col:
                raise SystemExit(f'ABORT: {path} has no {need!r} column; found {head[:8]}')
        cur_pos, block, prev = None, defaultdict(list), -1
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if f[col['#CHROM']] != want:
                continue
            pos = int(f[col['POS']])
            if pos != cur_pos:
                if cur_pos is not None:
                    yield cur_pos, block
                if cur_pos is not None and pos < prev:
                    raise SystemExit(f'ABORT: {path} is not position-sorted ({pos} < {prev})')
                prev, cur_pos, block = pos, pos, defaultdict(list)
            block[(f[col['REF']], f[col['ALT']])].append((
                f[col['effect']] or '.', f[col['impact']] or '.',
                f[col['gene']] or '.', f[col['geneid']] or '.', f[col['biotype']] or '.'))
        if cur_pos is not None:
            yield cur_pos, block


def collapse(rows, mode):
    """rows -> list of (effect, impact, gene, geneid, biotype), most-severe per unit."""
    if mode == 'variant_most_severe':
        best = min(rows, key=lambda r: (IMPACT_RANK.get(r[1], 99), r[0]))
        return [best]
    per = {}
    for eff, imp, gene, gid, bio in rows:
        k = gene
        r = IMPACT_RANK.get(imp, 99)
        if k not in per or (r, eff) < (IMPACT_RANK.get(per[k][1], 99), per[k][0]):
            per[k] = (eff, imp, gene, gid, bio)
    return sorted(per.values(), key=lambda r: (IMPACT_RANK.get(r[1], 99), r[2]))


def main():
    a = parse_args()
    eligible = {i.upper() for i in a.impacts}
    if 'MODIFIER' in eligible:
        raise SystemExit('ABORT: MODIFIER is never an eligible impact for a gene-based test.')

    keep = None
    if a.keep_variants:
        with open(a.keep_variants) as fh:
            keep = {ln.strip() for ln in fh if ln.strip()}

    stem = f'chr{a.chrom}'
    out_map = a.out_map or f'map.{stem}.tsv.gz'
    out_anno = a.out_anno or f'anno_counts.{stem}.tsv'
    out_qc = a.out_qc or f'map_qc.{stem}.tsv'

    idx = index_blocks(a.snpeff_index, a.chrom)
    idx_pos, idx_block = next(idx, (None, None))

    anno = Counter()          # (impact, effect, biotype) -> n variants, most-severe row only
    n_var = n_hit = n_miss = n_mod = n_kept = 0
    multi = 0

    with gzip.open(out_map, 'wt') as out:
        out.write('variant_id\tchrom\tpos\tref\talt\timpact\teffect\tbiotype\t'
                  'gene_symbol\tgene_id\tn_genes\n')
        for pos, vid, ref, alt in read_bim(a.bim, a.chrom, keep):
            n_var += 1
            while idx_pos is not None and idx_pos < pos:
                idx_pos, idx_block = next(idx, (None, None))
            rows = idx_block.get((ref, alt)) if (idx_pos == pos and idx_block) else None
            if not rows:
                n_miss += 1
                anno[(UNANNOTATED, '.', '.')] += 1
                continue
            n_hit += 1
            picked = collapse(rows, a.gene_assignment)
            # The distribution figure counts each variant ONCE, by its most
            # severe annotation, so the panels sum to the callset size.
            top = picked[0]
            anno[(top[1], top[0], top[4])] += 1
            usable = [r for r in picked if r[1] in eligible and r[2] not in ('.', '')]
            if not usable:
                n_mod += 1
                continue
            if len(usable) > 1:
                multi += 1
            for eff, imp, gene, gid, bio in usable:
                out.write(f'{vid}\t{a.chrom}\t{pos}\t{ref}\t{alt}\t{imp}\t{eff}\t{bio}\t'
                          f'{gene}\t{gid}\t{len(usable)}\n')
                n_kept += 1

    with open(out_anno, 'w') as fh:
        fh.write('Impact\tEffect\tBiotype\tCount\n')
        for (imp, eff, bio), n in sorted(anno.items(), key=lambda kv: -kv[1]):
            fh.write(f'{imp}\t{eff}\t{bio}\t{n}\n')

    with open(out_qc, 'w') as fh:
        fh.write('chrom\tn_variants_in\tn_annotated\tn_unannotated\tn_modifier_only\t'
                 'n_map_rows\tn_multigene_variants\tgene_assignment\timpacts\n')
        fh.write(f'{a.chrom}\t{n_var}\t{n_hit}\t{n_miss}\t{n_mod}\t{n_kept}\t{multi}\t'
                 f'{a.gene_assignment}\t{",".join(sorted(eligible))}\n')

    print(f'[build_gene_map] chr{a.chrom}: {n_var:,} variants in -> {n_hit:,} annotated '
          f'({n_miss:,} unannotated), {n_mod:,} MODIFIER-only, {n_kept:,} map rows '
          f'over {n_kept - multi if multi else n_kept:,}+ genes, {multi:,} multi-gene')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in build_gene_map: {e}', file=sys.stderr)
        sys.exit(1)
