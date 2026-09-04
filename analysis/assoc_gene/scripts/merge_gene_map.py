#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Merge the per-chromosome gene maps into the per-stratum artefacts
#           that FIX THE MULTIPLE-TESTING DENOMINATOR before any association is
#           run.
#
#           THE BARRIER. n_genes_mapped is decided here, once, from the map
#           alone. The >= --min-num-var filter is applied FIRST and the count of
#           surviving sets is what alpha is divided by, so the threshold cannot
#           drift with whatever an engine happened to return. Nothing downstream
#           re-derives it: map.<stratum>.tsv.gz is the FILTERED map, dropped sets
#           are absent from it and their variants are absent from
#           variants.<stratum>.txt, so every set handed to an engine is in the
#           denominator and every set in the denominator is handed to an engine.
#
#           SET NAMES AND THE CROSS-CHROMOSOME COLLISION. The canonical unit is
#           (gene_symbol, chrom), not the symbol. Small-RNA families (Y_RNA and
#           friends) really do reuse one symbol on several chromosomes in this
#           callset, and rvtest does not notice: it keys its set table on the
#           set name (Main.cpp:172, an OrderedMap) and silently MERGES the two
#           into one cross-chromosome set. A symbol living on more
#           than one chromosome therefore gets <symbol>@chr<N> on EVERY one of
#           its occurrences -- never a bare name for a colliding symbol -- so a
#           set name cannot mean one thing in one run and another in the next.
#           The split is decided from the union of all input map rows, not from
#           the stratum's subset, so the same gene carries the same name in every
#           stratum and a cross-stratum join stays honest.
#
#           POSITION COLLISIONS ARE ENUMERATED, NOT FIXED. rvtest expresses set
#           membership as 1-bp ranges, so a range at a position shared by two
#           records pulls in every record there; 0.52% of callset variants share
#           a position with another. This is a property of the range interface,
#           not a bug to repair here. map_collisions.<stratum>.tsv lists every
#           such position so SCAN can later classify an rvtest variant-count
#           mismatch as explained rather than unexplained.
# Component: assoc_gene
# ---------------------------------------------------------------------------
import argparse
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

# The pinned map schema. Read by name, never by position, but every input must
# carry all of it.
MAP_COLUMNS = ('variant_id', 'chrom', 'pos', 'ref', 'alt', 'impact', 'effect',
               'biotype', 'gene_symbol', 'gene_id', 'n_genes')
ANNO_COLUMNS = ('Impact', 'Effect', 'Biotype', 'Count')
VALID_IMPACTS = ('HIGH', 'MODERATE', 'LOW')
IMPACT_SLOT = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--map', nargs='+', action='append', required=True,
                   help='per-chromosome map.chr<N>.tsv.gz (repeatable)')
    p.add_argument('--anno', nargs='+', action='append', required=True,
                   help='per-chromosome anno_counts.chr<N>.tsv (repeatable)')
    p.add_argument('--qc', nargs='+', action='append', required=True,
                   help='per-chromosome map_qc.chr<N>.tsv (repeatable)')
    p.add_argument('--stratum', action='append', required=True, metavar='NAME=IMPACTS',
                   help='e.g. moderate_high=HIGH,MODERATE (repeatable)')
    p.add_argument('--min-num-var', type=int, default=2,
                   help='sets with fewer mapped variants are dropped BEFORE the denominator')
    p.add_argument('--alpha', type=float, default=0.05)
    p.add_argument('--collision-token', default='@',
                   help='separator for <symbol><token>chr<N> on multi-chromosome symbols')
    p.add_argument('--max-gene-span', type=int, default=3000000,
                   help='sets wider than this are flagged in map_qc.tsv')
    p.add_argument('--cohort', required=True)
    p.add_argument('--out-dir', default='.')
    return p.parse_args()


def flatten(nested):
    return [x for group in (nested or []) for x in group]


def opener(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def chrom_key(c):
    """Genomic order with the non-numeric contigs (X, Y, MT) after the autosomes."""
    return (0, int(c), '') if c.isdigit() else (1, 0, c)


def parse_strata(specs):
    out, seen = [], set()
    for spec in specs:
        name, sep, rest = spec.partition('=')
        name = name.strip()
        if not sep or not name:
            raise SystemExit(f'ABORT: --stratum {spec!r} is not NAME=IMPACTS')
        impacts = tuple(i.strip().upper() for i in rest.split(',') if i.strip())
        if not impacts:
            raise SystemExit(f'ABORT: --stratum {spec!r} names no impacts')
        if name in seen:
            raise SystemExit(f'ABORT: --stratum {name!r} given twice; the outputs would '
                             f'overwrite each other')
        for imp in impacts:
            if imp == 'MODIFIER':
                raise SystemExit('ABORT: MODIFIER is never a tested stratum -- it is what '
                                 'snpEff gives intergenic and deep-intronic variants, which a '
                                 'gene-based test has no principled way to assign.')
            if imp not in VALID_IMPACTS:
                raise SystemExit(f'ABORT: --stratum {name!r} names impact {imp!r}; expected '
                                 f'{"/".join(VALID_IMPACTS)}')
        seen.add(name)
        out.append((name, frozenset(impacts)))
    return out


def check_map_header(path, head):
    if 'set_name' in head:
        raise SystemExit(f'ABORT: {path} already carries a set_name column. That is a merged '
                         f'stratum map, not a per-chromosome map; re-running the merge on its '
                         f'own output would re-suffix the split sets.')
    col = {n: i for i, n in enumerate(head)}
    missing = [c for c in MAP_COLUMNS if c not in col]
    if missing:
        raise SystemExit(f'ABORT: {path} has no {missing} column(s); found {head}')
    return col


def read_header(path):
    with opener(path) as fh:
        return fh.readline().rstrip('\n').split('\t')


def iter_map(path, token):
    """Yield (raw_line, variant_id, chrom, pos, impact, gene_symbol) for one map file."""
    with opener(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        col = check_map_header(path, head)
        i_vid, i_chr = col['variant_id'], col['chrom']
        i_pos, i_imp, i_sym = col['pos'], col['impact'], col['gene_symbol']
        for ln, raw in enumerate(fh, 2):
            line = raw.rstrip('\n')
            if not line:
                continue
            f = line.split('\t')
            if len(f) != len(head):
                raise SystemExit(f'ABORT: {path}:{ln} has {len(f)} fields, expected {len(head)}')
            sym = f[i_sym]
            if token in sym:
                # A symbol carrying the token would make <symbol>@chr<N> ambiguous:
                # the merged set could no longer be split back into (symbol, chrom).
                raise SystemExit(
                    f'ABORT: {path}:{ln} gene_symbol {sym!r} already contains the collision '
                    f'token {token!r}. Re-run with --collision-token set to a character no '
                    f'snpEff symbol uses (and re-run every stratum, so no set name changes '
                    f'meaning between runs).')
            chrom = f[i_chr]
            if chrom.startswith('chr'):
                raise SystemExit(f'ABORT: {path}:{ln} chrom is {chrom!r}; the map carries the '
                                 f'PLAIN number, matching bim column 1')
            imp = f[i_imp]
            if imp not in VALID_IMPACTS:
                raise SystemExit(f'ABORT: {path}:{ln} impact {imp!r} is not one of '
                                 f'{"/".join(VALID_IMPACTS)}; the map must never carry MODIFIER')
            yield line, f[i_vid], chrom, int(f[i_pos]), imp, sym


def scan_maps(map_files, strata, token):
    """One pass: symbol -> chromosomes, and per-stratum (symbol, chrom) aggregates."""
    ref_head, chrom_owner, file_chroms = None, {}, {}
    sym_chroms = defaultdict(set)
    agg = {name: {} for name, _ in strata}
    n_rows = 0
    for path in map_files:
        head = read_header(path)
        check_map_header(path, head)
        if ref_head is None:
            ref_head = head
        elif head != ref_head:
            raise SystemExit(f'ABORT: {path} header differs from the first map file; the '
                             f'chromosomes were not all built by the same build_gene_map run')
        seen = set()
        for _line, _vid, chrom, pos, imp, sym in iter_map(path, token):
            n_rows += 1
            seen.add(chrom)
            sym_chroms[sym].add(chrom)
            for name, impacts in strata:
                if imp not in impacts:
                    continue
                k = (sym, chrom)
                e = agg[name].get(k)
                if e is None:
                    # build_gene_map keys its per-gene collapse on the symbol, so a
                    # (variant, symbol) pair is written at most once: the row count IS
                    # the distinct variant count and needs no dedup set.
                    agg[name][k] = e = [0, pos, pos, 0, 0, 0]
                e[0] += 1
                if pos < e[1]:
                    e[1] = pos
                if pos > e[2]:
                    e[2] = pos
                e[3 + IMPACT_SLOT[imp]] += 1
        for c in seen:
            if c in chrom_owner:
                raise SystemExit(f'ABORT: chromosome {c} appears in both {chrom_owner[c]} and '
                                 f'{path}; every variant would be counted twice')
            chrom_owner[c] = path
        file_chroms[path] = seen
    return ref_head, sym_chroms, agg, file_chroms, n_rows


def write_stratum(name, impacts, agg_sets, names, split_syms, order, header, args,
                  out_dir, token):
    """Filter, then write the five stratum artefacts. Returns the denominator row dict."""
    keep = {k: v for k, v in agg_sets.items() if v[0] >= args.min_num_var}
    dropped = len(agg_sets) - len(keep)

    idx_path = out_dir / f'gene_index.{name}.tsv'
    with open(idx_path, 'w') as fh:
        fh.write('set_name\tgene_symbol\tchrom\tn_var_map\tpos_min\tpos_max\tspan_bp\t'
                 'n_high\tn_moderate\tn_low\n')
        for (sym, chrom), v in sorted(keep.items(),
                                      key=lambda kv: (chrom_key(kv[0][1]), kv[1][1], kv[0][0])):
            fh.write(f'{names[(sym, chrom)]}\t{sym}\t{chrom}\t{v[0]}\t{v[1]}\t{v[2]}\t'
                     f'{v[2] - v[1]}\t{v[3]}\t{v[4]}\t{v[5]}\n')

    map_path = out_dir / f'map.{name}.tsv.gz'
    col_path = out_dir / f'map_collisions.{name}.tsv'
    vcount, vorder, n_colliding = {}, [], 0
    with gzip.open(map_path, 'wt') as gz, open(col_path, 'w') as cf:
        gz.write('\t'.join(header) + '\tset_name\n')
        cf.write('chrom\tpos\tn_set_names\tset_names\tvariant_ids\n')
        for path in order:
            # Held per input file, not per genome: bounded memory, and the flush order
            # follows the chromosome order of `order`.
            per_pos = {}
            for line, vid, chrom, pos, imp, sym in iter_map(path, token):
                if imp not in impacts:
                    continue
                k = (sym, chrom)
                if k not in keep:
                    continue
                sn = names[k]
                gz.write(f'{line}\t{sn}\n')
                seen = vcount.get(vid)
                if seen is None:
                    vcount[vid] = 1
                    vorder.append((chrom_key(chrom), pos, vid))
                else:
                    vcount[vid] = seen + 1
                e = per_pos.get((chrom, pos))
                if e is None:
                    per_pos[(chrom, pos)] = e = ({}, {})
                e[0][sn] = None      # dicts keep insertion order: stable output, no sets
                e[1][vid] = None
            for (chrom, pos) in sorted(per_pos, key=lambda cp: (chrom_key(cp[0]), cp[1])):
                sns, vids = per_pos[(chrom, pos)]
                if len(sns) < 2:
                    continue
                n_colliding += 1
                cf.write(f'{chrom}\t{pos}\t{len(sns)}\t{",".join(sorted(sns))}\t'
                         f'{",".join(sorted(vids))}\n')

    # Genomic order, not lexicographic: the list is read against a position-sorted bim.
    vorder.sort()
    with open(out_dir / f'variants.{name}.txt', 'w') as fh:
        for _k, _p, vid in vorder:
            fh.write(vid + '\n')

    n_mapped = len(keep)
    # 12 significant digits, not 6: gene_significance.py recomputes alpha/n_genes_mapped
    # and scan.py compares p-values against it, so the two must not disagree by a
    # rounding step at the threshold.
    thr = f'{args.alpha / n_mapped:.12g}' if n_mapped else 'NA'
    if not n_mapped:
        print(f'[merge_gene_map] WARNING: stratum {name!r} retained no set at '
              f'--min-num-var {args.min_num_var}; bonferroni_threshold is NA',
              file=sys.stderr)
    row = {
        'cohort': args.cohort, 'stratum': name, 'min_num_var': args.min_num_var,
        'n_genes_mapped': n_mapped, 'n_genes_dropped_lt_minvar': dropped,
        'n_variants': len(vcount),
        'n_variants_multigene': sum(1 for c in vcount.values() if c > 1),
        'n_sets_chrom_split': sum(1 for (s, _c) in keep if s in split_syms),
        'n_positions_colliding': n_colliding,
        'alpha': f'{args.alpha:g}', 'bonferroni_threshold': thr,
    }
    with open(out_dir / f'denominator.{name}.tsv', 'w') as fh:
        fh.write('\t'.join(row) + '\n')
        fh.write('\t'.join(str(v) for v in row.values()) + '\n')
    return keep, row


def merge_anno(anno_files, out_path):
    counts = Counter()
    for path in anno_files:
        with opener(path) as fh:
            head = fh.readline().rstrip('\n').split('\t')
            if tuple(head) != ANNO_COLUMNS:
                raise SystemExit(f'ABORT: {path} header is {head}, expected {list(ANNO_COLUMNS)}')
            for ln, raw in enumerate(fh, 2):
                line = raw.rstrip('\n')
                if not line:
                    continue
                f = line.split('\t')
                if len(f) != 4:
                    raise SystemExit(f'ABORT: {path}:{ln} has {len(f)} fields, expected 4')
                counts[(f[0], f[1], f[2])] += int(f[3])
    with open(out_path, 'w') as fh:
        fh.write('\t'.join(ANNO_COLUMNS) + '\n')
        for k in sorted(counts, key=lambda k: (-counts[k], k)):
            fh.write('\t'.join(k) + f'\t{counts[k]}\n')
    return sum(counts.values())


def read_qc(qc_files):
    """Read and cross-validate every per-chromosome QC row BEFORE anything is written."""
    ref_head, rows = None, []
    for path in qc_files:
        with opener(path) as fh:
            head = fh.readline().rstrip('\n').split('\t')
            if ref_head is None:
                ref_head = head
            elif head != ref_head:
                raise SystemExit(f'ABORT: {path} header differs from the first --qc file')
            body = [ln.rstrip('\n').split('\t') for ln in fh if ln.strip()]
        if not body:
            raise SystemExit(f'ABORT: {path} has no data row')
        rows.extend(body)
    if 'chrom' not in ref_head:
        raise SystemExit(f'ABORT: --qc files have no chrom column; found {ref_head}')
    i_chrom = ref_head.index('chrom')
    for name in ('gene_assignment', 'impacts'):
        if name in ref_head:
            i = ref_head.index(name)
            vals = sorted({r[i] for r in rows})
            if len(vals) > 1:
                # Mixing chromosomes built under different rules would silently mix
                # two definitions of "the map" inside one denominator.
                raise SystemExit(f'ABORT: --qc files disagree on {name}: {vals}. The '
                                 f'chromosomes were not built by one build_gene_map policy.')
    rows.sort(key=lambda r: chrom_key(r[i_chrom]))
    return ref_head, rows


def write_qc(ref_head, rows, out_path, over, max_span):
    i_chrom = ref_head.index('chrom')
    total = []
    for i, cname in enumerate(ref_head):
        if i == i_chrom:
            total.append('TOTAL')
            continue
        vals = [r[i] for r in rows]
        try:
            total.append(str(sum(int(v) for v in vals)))
        except ValueError:
            u = sorted(set(vals))
            total.append(u[0] if len(u) == 1 else 'MIXED')

    with open(out_path, 'w') as fh:
        fh.write('\t'.join(ref_head) + '\tmax_gene_span\tn_sets_over_max_span\t'
                 'sets_over_max_span\n')
        n_over_total = 0
        for r in rows:
            hits = over.get(r[i_chrom], {})
            n_over_total += len(hits)
            detail = ';'.join(
                # '+' joins the strata: '|' is a plausible --collision-token and
                # would then be ambiguous inside a set name.
                f'{sn}:{span}:{"+".join(sorted(strata))}'
                for sn, (span, strata) in sorted(hits.items(), key=lambda kv: -kv[1][0])) or '.'
            fh.write('\t'.join(r) + f'\t{max_span}\t{len(hits)}\t{detail}\n')
        fh.write('\t'.join(total) + f'\t{max_span}\t{n_over_total}\t.\n')
    return n_over_total


def main():
    args = parse_args()
    strata = parse_strata(args.stratum)
    map_files, anno_files, qc_files = flatten(args.map), flatten(args.anno), flatten(args.qc)
    if not (len(map_files) == len(anno_files) == len(qc_files)):
        raise SystemExit(f'ABORT: {len(map_files)} --map, {len(anno_files)} --anno and '
                         f'{len(qc_files)} --qc file(s); one of each per chromosome is '
                         f'required, otherwise a chromosome is silently absent from a summary')
    if args.min_num_var < 1:
        raise SystemExit(f'ABORT: --min-num-var {args.min_num_var} < 1')
    if not 0 < args.alpha < 1:
        raise SystemExit(f'ABORT: --alpha {args.alpha} is not in (0, 1)')
    if not args.collision_token:
        raise SystemExit('ABORT: --collision-token must be a non-empty string')
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    token = args.collision_token

    # Everything that can abort reads before anything is written: a run that dies on
    # a bad --qc file must not leave half a stratum's artefacts on disk for a
    # downstream step to pick up.
    qc_head, qc_rows = read_qc(qc_files)
    head, sym_chroms, agg, file_chroms, n_rows = scan_maps(map_files, strata, token)

    # Decided over the union of all input rows, so the name of a gene does not
    # depend on which strata were requested.
    split_syms = {s for s, cs in sym_chroms.items() if len(cs) > 1}
    names = {}
    for sym, chroms in sym_chroms.items():
        for chrom in chroms:
            names[(sym, chrom)] = f'{sym}{token}chr{chrom}' if sym in split_syms else sym
    order = sorted(map_files,
                   key=lambda p: min((chrom_key(c) for c in file_chroms[p]), default=(2, 0, '')))

    print(f'[merge_gene_map] {len(map_files)} chromosome(s), {n_rows:,} map row(s), '
          f'{len(sym_chroms):,} symbol(s) -> {len(names):,} (symbol, chrom) set(s)')
    if split_syms:
        shown = ','.join(sorted(split_syms)[:8]) + ('...' if len(split_syms) > 8 else '')
        print(f'    {len(split_syms):,} symbol(s) on >1 chromosome, suffixed '
              f'{token}chr<N> everywhere: {shown}')

    over = defaultdict(dict)   # chrom -> set_name -> [max span over strata, {strata}]
    dens = []
    for name, impacts in strata:
        keep, row = write_stratum(name, impacts, agg[name], names, split_syms,
                                  order, head, args, out_dir, token)
        dens.append(row)
        for (sym, chrom), v in keep.items():
            span = v[2] - v[1]
            if span > args.max_gene_span:
                e = over[chrom].setdefault(names[(sym, chrom)], [span, set()])
                e[0] = max(e[0], span)
                e[1].add(name)
        print(f'    {name:>20s}: {row["n_genes_mapped"]:,} set(s) tested '
              f'(+{row["n_genes_dropped_lt_minvar"]:,} dropped at <{args.min_num_var} var), '
              f'{row["n_variants"]:,} variant(s), '
              f'{row["n_sets_chrom_split"]:,} chrom-split, '
              f'{row["n_positions_colliding"]:,} colliding position(s), '
              f'p < {row["bonferroni_threshold"]}')
        top = sorted(keep.items(), key=lambda kv: (-kv[1][0], -(kv[1][2] - kv[1][1])))[:10]
        for (sym, chrom), v in top:
            flag = ' OVER-SPAN' if (v[2] - v[1]) > args.max_gene_span else ''
            print(f'        {names[(sym, chrom)]:<24s} n_var={v[0]:<6d} '
                  f'span={v[2] - v[1]:,} bp{flag}')

    n_over = write_qc(qc_head, qc_rows, out_dir / 'map_qc.tsv', over, args.max_gene_span)
    n_anno = merge_anno(anno_files, out_dir / 'anno_counts.tsv')
    print(f'[merge_gene_map] anno_counts.tsv: {n_anno:,} variant(s) counted once each; '
          f'map_qc.tsv: {len(qc_rows)} chromosome row(s) + TOTAL, '
          f'{n_over:,} set(s) over --max-gene-span {args.max_gene_span:,}')
    print(f'[merge_gene_map] denominator fixed for {len(dens)} stratum/strata -> {out_dir}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in merge_gene_map: {e}', file=sys.stderr)
        sys.exit(1)
