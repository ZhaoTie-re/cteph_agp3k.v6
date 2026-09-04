#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Concatenate the 22 per-chromosome rvtest tables for one
#           (cohort, stratum, method) into one table on the gene_scan schema,
#           joined to the per-gene counts, with the effect size derived once.
#
#           THE COLUMN THAT IS NOT THE COLUMN. The pipeline drives rvtest with
#           --setFile, not --geneFile, and that renames the output's first
#           column from `Gene` to `Range` -- which holds the SET NAME -- while
#           adding a second column `RANGE`, the consolidated 1-bp range list,
#           ~30 kB of text for a large gene. The two differ only in case. Taken
#           by position, or by name through any case-folding path, 30 kB of
#           coordinates lands where the set name belongs, joins to nothing, and
#           the whole scan comes back unmapped with no error anywhere. So the
#           header is ASSERTED, `Range` becomes set_name at the point of read,
#           and `RANGE` is dropped before one row is kept.
#
#           WHY THE ASSERT IS FATAL. The one way to see a `Gene` first column
#           here is for someone to have put --geneFile back, and then the tested
#           sets are rvtest's own refFlat transcript spans -- not the snpEff
#           map, and not what MERGE_MAP counted the denominator from.
#
#           THE EFFECT SIZE. rvtest's --burden cmc is a SCORE test: a p-value
#           and no coefficient. The OR comes from --burden cmcWald, run as a
#           second model in the SAME invocation on the SAME collapsed genotype
#           and covariates, so OR = exp(beta) is that model's own estimate. Two
#           consequences are stated rather than left to be discovered: the Wald
#           p is NOT published (two p-values for one model on one row invites
#           the reader to pick one), and 2*Phi(-|beta/se|) will not reproduce
#           the published score p. The CI is the Wald normal approximation,
#           exp(beta -/+ 1.96 se). Kernel rows (skato) carry no effect size:
#           SKAT-O does not estimate one. rho, the mixing parameter rvtest
#           chose, is kept on skato rows.
#
#           WHY THE PARTS ARE NOT TRUSTED TO BE IN ORDER. Nextflow stages them
#           as part_1..part_22 in staging order, which is not chromosome order.
#           Rows are sorted into genomic order here from the joined index, so
#           merged.tsv is byte-identical whatever order the parts arrived in.
# Component: assoc_gene
# ---------------------------------------------------------------------------
import argparse
import math
import sys

# The pinned rvtest --setFile header. The first five are asserted exactly; the
# tail identifies which test actually ran and is checked against --method, so a
# cmc table cannot be merged under the skato label.
SETFILE_PREFIX = ('Range', 'RANGE', 'N_INFORMATIVE', 'NumVar', 'NumPolyVar')
METHOD_TAIL = {'cmc': ('NonRefSite', 'Pvalue'), 'skato': ('Q', 'rho', 'Pvalue')}
# The second model of the same rvtest invocation, present only for --method cmc.
# Measured header (rvtests 20190205, --burden cmc,cmcWald --setFile):
#   Range RANGE N_INFORMATIVE NumVar NumPolyVar NonRefSite Beta SE Pvalue
WALD_TAIL = ('NonRefSite', 'Beta', 'SE', 'Pvalue')

INDEX_COLUMNS = ('set_name', 'gene_symbol', 'chrom', 'n_var_map', 'pos_min', 'pos_max',
                 'span_bp', 'n_high', 'n_moderate', 'n_low')
SET_STATS_READ = ('set_name', 'mac', 'mac_case', 'mac_control',
                  'n_carrier_case', 'n_carrier_control')

# gene_scan.tsv, in order. The last four columns stay NA here: both decision
# rules are properties of the DENOMINATOR, which this step is deliberately not
# handed. SCAN owns the denominator and fills all four.
SCAN_COLUMNS = ('cohort', 'stratum', 'method', 'set_name', 'gene_symbol', 'chrom',
                'pos_min', 'pos_max', 'n_var_map', 'n_var_engine', 'n_informative',
                'pvalue', 'rho', 'beta', 'se', 'or', 'or_l95', 'or_u95',
                'mac', 'mac_case', 'mac_control', 'n_carrier_case', 'n_carrier_control',
                'bonferroni_threshold', 'significant_bonferroni', 'fdr_bh', 'significant_bh')

Z_95 = 1.959963984540054
NA = 'NA'
MAX_REPORT = 12


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assoc', nargs='+', required=True,
                   help='per-chromosome rvtest chr<N>.assoc, staged as part_*.assoc')
    p.add_argument('--wald', nargs='+', required=True,
                   help='the matching chr<N>.wald.assoc parts. For --method cmc these are '
                        'rvtest --burden cmcWald output and supply beta/se; for a kernel '
                        'method they must all be EMPTY, which is asserted rather than '
                        'assumed. One part per --assoc part.')
    p.add_argument('--n-covar', type=int, required=True,
                   help='number of --covar-name columns the RVTEST task passed. cmcWald '
                        'writes 1 + n_covar rows per gene (Model.h:961 loops over '
                        'X = [Intercept | collapsed | covariates]) and only the FIRST is '
                        'the burden term; passing the count ties the .nf covariate list to '
                        'this assumption instead of inferring it from the file alone.')
    p.add_argument('--gene-index', required=True, help='gene_index.<stratum>.tsv')
    p.add_argument('--set-stats', required=True,
                   help='set_stats.tsv from EMIT_TEST_INPUTS; per-set allele and carrier counts')
    p.add_argument('--cohort', required=True)
    p.add_argument('--stratum', required=True)
    p.add_argument('--method', required=True, choices=sorted(METHOD_TAIL))
    p.add_argument('--out', default='merged.tsv')
    return p.parse_args()


def chrom_key(c):
    """Genomic order with the non-numeric contigs after the autosomes, NA last."""
    s = str(c)
    if s.isdigit():
        return (0, int(s), '')
    if s in ('', NA):
        return (2, 0, '')
    return (1, 0, s)


def warn(msg):
    """stderr, but only after stdout has caught up -- the two are buffered
    differently once Nextflow captures them and the log otherwise reads out of order."""
    sys.stdout.flush()
    print(msg, file=sys.stderr, flush=True)


def sample(items):
    """A short, deterministic excerpt for a report line."""
    s = sorted(items)
    more = f' ... (+{len(s) - MAX_REPORT} more)' if len(s) > MAX_REPORT else ''
    return ', '.join(s[:MAX_REPORT]) + more


def check_header(path, head, method):
    """Fatal unless this is the --setFile header for the method we were told ran."""
    if tuple(head[:5]) != SETFILE_PREFIX:
        raise SystemExit(
            f'ABORT: {path} starts with {head[:5]}, expected {list(SETFILE_PREFIX)}.\n'
            f'       A first column named "Gene" means rvtest ran under --geneFile, so the '
            f'tested sets are its own refFlat transcript spans -- not the snpEff map, and '
            f'not what the denominator was counted from. Nothing downstream can detect '
            f'that from the table alone, so it is refused here.')
    tail, expected = tuple(head[5:]), METHOD_TAIL[method]
    if tail == expected:
        return
    other = [m for m, t in METHOD_TAIL.items() if t == tail]
    if other:
        raise SystemExit(f'ABORT: {path} carries the {other[0]} tail {list(tail)} but was '
                         f'merged as {method}. The p-values would be labelled with a test '
                         f'that did not produce them.')
    raise SystemExit(f'ABORT: {path} tail is {list(tail)}, expected {list(expected)} for '
                     f'{method}; this is not the rvtest output this merge was written for.')


def read_gene_index(path):
    """set_name -> {gene_symbol, chrom, pos_min, pos_max, n_var_map}."""
    idx = {}
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        col = {n: i for i, n in enumerate(head)}
        missing = [c for c in INDEX_COLUMNS if c not in col]
        if missing:
            raise SystemExit(f'ABORT: {path} has no {missing} column(s); found {head}')
        for ln, raw in enumerate(fh, 2):
            line = raw.rstrip('\n')
            if not line:
                continue
            f = line.split('\t')
            if len(f) != len(head):
                raise SystemExit(f'ABORT: {path}:{ln} has {len(f)} fields, expected {len(head)}')
            sn = f[col['set_name']]
            if sn in idx:
                raise SystemExit(f'ABORT: {path}:{ln} lists set {sn!r} twice; the index is the '
                                 f'denominator, so a set counted twice inflates it')
            idx[sn] = {k: f[col[k]] for k in
                       ('gene_symbol', 'chrom', 'pos_min', 'pos_max', 'n_var_map')}
    if not idx:
        raise SystemExit(f'ABORT: {path} lists no set; there is nothing to join to')
    return idx


def read_assoc(path, method):
    """[(set_name, n_informative, n_var_engine, pvalue, rho)], plus the contig the ranges name."""
    rows, contigs = [], set()
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        check_header(path, head, method)
        col = {n: i for i, n in enumerate(head)}
        i_set, i_rng = col['Range'], col['RANGE']
        i_inf, i_num, i_p = col['N_INFORMATIVE'], col['NumVar'], col['Pvalue']
        i_rho = col.get('rho')
        for ln, raw in enumerate(fh, 2):
            line = raw.rstrip('\n')
            if not line:
                continue
            f = line.split('\t')
            if len(f) != len(head):
                raise SystemExit(f'ABORT: {path}:{ln} has {len(f)} fields, expected {len(head)}')
            # The only read of RANGE in this component: the contig prefix of its
            # first range, sliced rather than split so the ~30 kB cell is never
            # copied. Everything after it is discarded with the field.
            cell = f[i_rng]
            cut = cell.find(':')
            if cut > 0:
                contigs.add(cell[:cut].removeprefix('chr'))
            rows.append((f[i_set], f[i_inf], f[i_num], f[i_p],
                         f[i_rho] if i_rho is not None else NA))
    return rows, contigs


def read_wald(path, n_covar):
    """set_name -> (beta, se) from a --burden cmcWald part, or {} if the part is empty.

    cmcWald writes 1 + n_covar rows per gene, all carrying the SAME Range and the
    same ~30 kB RANGE cell: Model.h:961 loops `for (i = 1; i < X.cols; ++i)` over
    X = [Intercept | collapsed | covariates], so i = 1 is the collapsed burden and
    the rest are the covariate coefficients. Only column LABELS would distinguish
    them and rvtest does not write any, so the burden row is identified
    POSITIONALLY -- it is the first row of each block. Verified against a real run
    with two covariates: three rows per gene, the first agreeing with the CMC
    score p to 3 significant digits (0.25879 Wald vs 0.25805 score) and the other
    two being the covariates.

    The block height is asserted twice: against --n-covar, which ties this to the
    .nf covariate list, and against the file's own constant height, which catches
    a file that is internally inconsistent. Neither alone is sufficient.
    """
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        if not head or head == ['']:
            return {}
        if tuple(head[:5]) != SETFILE_PREFIX or tuple(head[5:]) != WALD_TAIL:
            raise SystemExit(f'ABORT: {path} header is {head}, expected '
                             f'{list(SETFILE_PREFIX) + list(WALD_TAIL)}. This is not a '
                             f'--burden cmcWald --setFile output.')
        col = {n: i for i, n in enumerate(head)}
        i_set, i_b, i_se = col['Range'], col['Beta'], col['SE']
        order, blocks = [], {}
        for ln, raw in enumerate(fh, 2):
            line = raw.rstrip('\n')
            if not line:
                continue
            f = line.split('\t')
            if len(f) != len(head):
                raise SystemExit(f'ABORT: {path}:{ln} has {len(f)} field(s), expected '
                                 f'{len(head)}')
            sn = f[i_set]
            if sn not in blocks:
                blocks[sn] = []
                order.append(sn)
            elif order[-1] != sn:
                raise SystemExit(f'ABORT: {path}:{ln} returns to set {sn!r} after '
                                 f'{order[-1]!r}; cmcWald blocks must be contiguous or the '
                                 f'burden row cannot be identified positionally.')
            blocks[sn].append((f[i_b], f[i_se]))
    if not blocks:
        raise SystemExit(f'ABORT: {path} has a cmcWald header and no data row')
    heights = {len(v) for v in blocks.values()}
    if len(heights) != 1:
        raise SystemExit(f'ABORT: {path} block heights are {sorted(heights)}; cmcWald writes '
                         f'1 + n_covar rows for EVERY gene.')
    height = heights.pop()
    if height != 1 + n_covar:
        raise SystemExit(
            f'ABORT: {path} has {height} row(s) per gene but --n-covar {n_covar} implies '
            f'{1 + n_covar} (one burden term plus one per covariate). Either the RVTEST '
            f'task passed a different --covar-name list than this merge was told, or rvtest '
            f'dropped a covariate. The burden row is taken POSITIONALLY, so publishing a '
            f'beta under a mismatched block height risks publishing a covariate\'s.')
    return {sn: rows[0] for sn, rows in blocks.items()}


def read_set_stats(path):
    """set_name -> the count columns, from EMIT_TEST_INPUTS."""
    out = {}
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        missing = [c for c in SET_STATS_READ if c not in head]
        if missing:
            raise SystemExit(f'ABORT: {path} has no {missing} column(s); found {head}. '
                             f'This is not set_stats.tsv.')
        col = {n: i for i, n in enumerate(head)}
        for ln, raw in enumerate(fh, 2):
            line = raw.rstrip('\n')
            if not line:
                continue
            f = line.split('\t')
            if len(f) != len(head):
                raise SystemExit(f'ABORT: {path}:{ln} has {len(f)} field(s), expected '
                                 f'{len(head)}')
            sn = f[col['set_name']]
            if sn in out:
                raise SystemExit(f'ABORT: {path} lists set {sn!r} twice')
            out[sn] = {k: f[col[k]] for k in SET_STATS_READ if k != 'set_name'}
    if not out:
        raise SystemExit(f'ABORT: {path} lists no set')
    return out


def as_int(v):
    """rvtest writes NA for a set it could not evaluate; keep that, never a 0."""
    try:
        return str(int(v))
    except (TypeError, ValueError):
        return NA


def as_pvalue(v):
    """Kept VERBATIM when it parses -- a float round-trip loses digits that matter at 1e-300."""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return NA, False
    if x != x or x in (float('inf'), float('-inf')) or not 0.0 <= x <= 1.0:
        return NA, False
    return v, True


def clean_float(v):
    """rvtest writes a plain float or NA (ModelFitter.h:46 resets to defaultValue
    "NA" when a fit fails), so a failed Wald fit gives NA here and never a stale
    value from the previous gene."""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if x != x or x in (float('inf'), float('-inf')) else x


def effect_size(beta_s, se_s):
    """(beta, se, or, or_l95, or_u95) as strings, all NA unless both parse and se > 0."""
    b, s = clean_float(beta_s), clean_float(se_s)
    if b is None or s is None or s <= 0:
        return NA, NA, NA, NA, NA
    try:
        lo, hi = math.exp(b - Z_95 * s), math.exp(b + Z_95 * s)
        o = math.exp(b)
    except OverflowError:
        return beta_s, se_s, NA, NA, NA
    return beta_s, se_s, f'{o:.6g}', f'{lo:.6g}', f'{hi:.6g}'


def main():
    args = parse_args()
    idx = read_gene_index(args.gene_index)
    stats = read_set_stats(args.set_stats)

    # ---- the second model of the same invocation, when this method has one.
    if len(args.wald) != len(args.assoc):
        raise SystemExit(f'ABORT: {len(args.wald)} --wald part(s) for {len(args.assoc)} '
                         f'--assoc part(s); RVTEST emits one of each per chromosome.')
    wald, n_empty = {}, 0
    for wp in args.wald:
        got = read_wald(wp, args.n_covar)
        if not got:
            n_empty += 1
            continue
        dupe = set(got) & set(wald)
        if dupe:
            raise SystemExit(f'ABORT: {wp} repeats set(s) {sample(dupe)} already seen in '
                             f'another wald part')
        wald.update(got)
    # A KERNEL METHOD ESTIMATES NO EFFECT SIZE. RVTEST writes an empty
    # chr<N>.wald.assoc for it deliberately (params.RvtestMethods wald: null),
    # so emptiness is ASSERTED per method rather than inferred.
    if args.method == 'cmc' and n_empty:
        raise SystemExit(f'ABORT: {n_empty} of {len(args.wald)} wald part(s) are empty for '
                         f'--method cmc, which runs --burden cmc,cmcWald.')
    if args.method != 'cmc' and n_empty != len(args.wald):
        raise SystemExit(f'ABORT: --method {args.method} estimates no effect size, so every '
                         f'wald part must be empty; {len(args.wald) - n_empty} carried data.')

    merged, seen, n_bad_p, n_no_or = [], {}, 0, 0
    per_chrom = {}
    for path in args.assoc:
        rows, contigs = read_assoc(path, args.method)
        if len(contigs) > 1:
            raise SystemExit(f'ABORT: {path} names contigs {sorted(contigs)}; one part is one '
                             f'chromosome. A set spanning two chromosomes is silently merged '
                             f'by rvtest.')
        file_chrom = contigs.pop() if contigs else NA
        for set_name, n_inf, n_var, pval, rho in rows:
            if set_name in seen:
                raise SystemExit(f'ABORT: set {set_name!r} is returned by both {seen[set_name]} '
                                 f'and {path}. gene_scan.tsv is one row per gene x method.')
            seen[set_name] = path
            g = idx.get(set_name)
            st = stats.get(set_name)
            p, ok = as_pvalue(pval)
            n_bad_p += not ok
            chrom = g['chrom'] if g else file_chrom
            if args.method == 'cmc':
                w = wald.get(set_name)
                beta, se, o, lo, hi = effect_size(*w) if w else (NA,) * 5
                n_no_or += (o == NA)
                rho_s = NA
            else:
                beta = se = o = lo = hi = NA
                rho_s = rho if clean_float(rho) is not None else NA
            merged.append({
                'cohort': args.cohort, 'stratum': args.stratum, 'method': args.method,
                'set_name': set_name,
                'gene_symbol': g['gene_symbol'] if g else NA,
                'chrom': chrom,
                'pos_min': g['pos_min'] if g else NA,
                'pos_max': g['pos_max'] if g else NA,
                'n_var_map': g['n_var_map'] if g else NA,
                'n_var_engine': as_int(n_var), 'n_informative': as_int(n_inf),
                'pvalue': p, 'rho': rho_s,
                'beta': beta, 'se': se, 'or': o, 'or_l95': lo, 'or_u95': hi,
                # The counts describe the VARIANT SET, not the test, so they go
                # on both method rows.
                'mac': st['mac'] if st else NA,
                'mac_case': st['mac_case'] if st else NA,
                'mac_control': st['mac_control'] if st else NA,
                'n_carrier_case': st['n_carrier_case'] if st else NA,
                'n_carrier_control': st['n_carrier_control'] if st else NA,
                'bonferroni_threshold': NA, 'significant_bonferroni': NA,
                'fdr_bh': NA, 'significant_bh': NA,
            })
            per_chrom[chrom] = per_chrom.get(chrom, 0) + 1
        print(f'[merge_rvtest] {path}: chr{file_chrom}, {len(rows):,} row(s)')

    if not merged:
        raise SystemExit(f'ABORT: {len(args.assoc)} part file(s) carried no data row. rvtest '
                         f'exits 0 on an empty result, so this is not visible upstream.')

    # ---- both models ran over the same sets, and so did the counts.
    returned = set(seen)
    if wald:
        only_score, only_wald = returned - set(wald), set(wald) - returned
        if only_score or only_wald:
            raise SystemExit(
                f'ABORT: the score and Wald models of the same invocation returned different '
                f'set(s): {len(only_score):,} in {args.method} only [{sample(only_score)}], '
                f'{len(only_wald):,} in cmcWald only [{sample(only_wald)}].')
    no_stats = returned - set(stats)
    if no_stats:
        raise SystemExit(
            f'ABORT: {len(no_stats):,} set(s) rvtest returned carry no row in '
            f'{args.set_stats} [{sample(no_stats)}]. The counts and the tested sets come '
            f'from the same EMIT_TEST_INPUTS task, so a gap means a stale set_stats.tsv '
            f'was staged against a fresh scan.')

    def sort_key(r):
        pos = int(r['pos_min']) if r['pos_min'].isdigit() else -1
        return (chrom_key(r['chrom']), pos, r['set_name'])

    merged.sort(key=sort_key)
    with open(args.out, 'w') as fh:
        fh.write('\t'.join(SCAN_COLUMNS) + '\n')
        for r in merged:
            fh.write('\t'.join(r[c] for c in SCAN_COLUMNS) + '\n')

    counts = ' '.join(f'chr{c}:{per_chrom[c]:,}' for c in sorted(per_chrom, key=chrom_key))
    print(f'[merge_rvtest] per-chromosome rows: {counts}')
    print(f'[merge_rvtest] {args.cohort}/{args.stratum}/{args.method}: {len(merged):,} set(s) '
          f'from {len(args.assoc)} part file(s) -> {args.out}')

    missing = set(idx) - set(seen)
    extra = set(seen) - set(idx)
    if missing:
        print(f'[merge_rvtest] {len(missing):,} of {len(idx):,} indexed set(s) absent from the '
              f'rvtest output: {sample(missing)}')
    if extra:
        warn(f'[merge_rvtest] WARNING: {len(extra):,} returned set(s) are absent from '
             f'{args.gene_index} and carry no gene_symbol/pos: {sample(extra)}. rvtest cannot '
             f'invent a set name -- the set file and the gene index are not from the same '
             f'EMIT_TEST_INPUTS task.')
    idx_chroms = {v['chrom'] for v in idx.values()}
    empty = sorted(idx_chroms - set(per_chrom), key=chrom_key)
    if empty:
        warn(f'[merge_rvtest] WARNING: chromosome(s) {", ".join(empty)} carry indexed sets but '
             f'returned no row at all; a whole rvtest part is missing or empty.')
    if n_bad_p:
        print(f'[merge_rvtest] {n_bad_p:,} row(s) had an unusable Pvalue (NA/nan/out of '
              f'[0,1]) and were written as NA')
    if args.method == 'cmc' and n_no_or:
        print(f'[merge_rvtest] {n_no_or:,} cmc row(s) have no OR: rvtest reset the Wald fit '
              f'to NA (typically a set with no minor allele in one group)')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in merge_rvtest: {e}', file=sys.stderr)
        sys.exit(1)
