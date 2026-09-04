#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : BARRIER B3. Both methods' tables for one (cohort, stratum) become
#           one scan: the map identity is proved per gene, and the two decision
#           rules are applied over the family the denominator fixed.
#
#           THE MAP IDENTITY IS  n_var_map == NumVar, OUTER-JOINED over the
#           union of the gene index and rvtest's output. An inner join would
#           make the one failure that matters -- a gene that was mapped and
#           never returned -- vanish from the very table written to detect it.
#           So every set gets a row in map_identity.tsv whether it passed,
#           failed, or was never returned.
#
#           WHY THE IDENTITY CAN HOLD EXACTLY EVEN WHEN NumVar > n_var_map.
#           rvtest is handed one 1-bp range per variant, and a tabix interval
#           matches a record over [POS, POS+len(REF)-1]. A multi-base-REF
#           record (a deletion) is therefore matched a SECOND time by the range
#           of a later mapped variant that sits inside its span, and rvtest
#           counts it twice. span_overlap_excess() predicts that count exactly
#           -- 81,746 / 81,746 genes across nine cohort x stratum cells on the
#           reference data -- so such a gene is classified
#           explained_by_span_overlap and the check stays an EQUALITY rather
#           than being loosened to an inequality. This is also a real analytic
#           fact, not only bookkeeping: rvtest genuinely enters the deletion
#           twice. Harmless for CMC (a binary collapse), but the SKAT-O kernel
#           double-weights it; the flag is carried so a reader can see which
#           genes it touches.
#
#           WHY THE THRESHOLD IS READ AND NEVER RECOMPUTED. alpha/n_genes_mapped
#           is fixed at barrier B1, before a single association ran. The STRING
#           from denominator.<stratum>.tsv is copied verbatim into every row,
#           and the denominator's own arithmetic is confirmed rather than
#           trusted.
#
#           WHY n_genes_missing IS SURFACED RATHER THAN SUBTRACTED. The
#           correction is for the experiment that was designed, not for the
#           subset that converged, so n_genes_mapped stays the denominator and
#           the shortfall gets its own column.
#
#           WHY lambda_gc IS ALSO REPORTED BY n_var. The statistic is a discrete
#           count over a handful of minor alleles, and for a two-variant gene
#           the normal approximation both tests rely on is anti-conservative
#           near p = 0.5 (a pile-up just under the median) while being fine in
#           the tail. That inflates the median-based lambda for the n_var = 2
#           genes and only them, so lambda is reported overall AND split at
#           n_var = 2 / >= 3. lambda ~ 1 is not the target here; a null
#           simulation on the real MAC distribution gives ~0.6 for an exact
#           test. docs/METHODS.md carries the measurement.
#
#           WHY lambda_gc GOES THROUGH inv_cdf(p/2) AND NOT inv_cdf(1 - p/2).
#           In binary64, 1 - p/2 is exactly 1.0 for every p below ~2.2e-16 and
#           the quantile there is +inf. p/2 is exact down to the subnormals.
# Component: assoc_gene
# ---------------------------------------------------------------------------
import argparse
import gzip
import statistics
import sys
from statistics import NormalDist

# gene_scan.tsv, pinned. MERGE_RVTEST writes this header, so it is asserted here
# EXACTLY rather than by name: a merged.tsv whose columns moved is a sibling that
# was edited without this barrier being told.
SCAN_COLUMNS = ('cohort', 'stratum', 'method', 'set_name', 'gene_symbol', 'chrom',
                'pos_min', 'pos_max', 'n_var_map', 'n_var_engine', 'n_informative',
                'pvalue', 'rho', 'beta', 'se', 'or', 'or_l95', 'or_u95',
                'mac', 'mac_case', 'mac_control', 'n_carrier_case', 'n_carrier_control',
                'bonferroni_threshold', 'significant_bonferroni', 'fdr_bh', 'significant_bh')

QC_COLUMNS = ('cohort', 'stratum', 'method', 'n_genes_mapped', 'n_genes_returned',
              'n_genes_missing', 'lambda_gc', 'lambda_gc_nvar2', 'lambda_gc_nvar_ge3',
              'min_p', 'n_significant_bonferroni', 'n_significant_bh', 'bh_threshold',
              'bonferroni_threshold', 'min_num_var')

IDENTITY_COLUMNS = ('cohort', 'stratum', 'chrom', 'set_name', 'n_var_map', 'n_var_engine',
                    'span_overlap_excess', 'ok', 'reason')

INDEX_COLUMNS = ('set_name', 'gene_symbol', 'chrom', 'n_var_map', 'pos_min', 'pos_max',
                 'span_bp', 'n_high', 'n_moderate', 'n_low')

DENOM_COLUMNS = ('cohort', 'stratum', 'min_num_var', 'n_genes_mapped',
                 'n_genes_dropped_lt_minvar', 'n_variants', 'n_variants_multigene',
                 'n_sets_chrom_split', 'n_positions_colliding', 'alpha',
                 'bonferroni_threshold')

COLLISION_COLUMNS = ('chrom', 'pos', 'n_set_names', 'set_names', 'variant_ids')

# The closed vocabulary. A typo'd method would not error anywhere: it would add a
# third scan_qc row that SIGNIFICANCE never asks for and remove one it does.
METHODS = ('cmc', 'skato')

REASONS = ('ok', 'explained_by_span_overlap', 'explained_by_collision', 'missing',
           'unexplained')
FAIL_REASONS = ('missing', 'unexplained')

# median(chi2_1) under the null. Pinned by the spec, not recomputed from
# inv_cdf(0.25)**2 -- that lands 2 ulp away and would make lambda_gc differ in
# the last digits between this script and gene_significance.py.
CHI2_MEDIAN_1DF = 0.4549364231195729
P_FLOOR = 1e-300
# The tolerance on alpha / n against the denominator file is on MERGE_MAP's
# 12-significant-digit round trip, the same 1e-9 gene_significance.py uses.
REL_TOL = 1e-9

NA = 'NA'
MAX_REPORT = 20


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assoc', nargs='+', required=True,
                   help='merged.tsv from MERGE_RVTEST, one per method')
    p.add_argument('--gene-index', required=True, help='gene_index.<stratum>.tsv')
    p.add_argument('--denominator', required=True,
                   help='denominator.<stratum>.tsv; the threshold is READ from it')
    p.add_argument('--collisions', required=True, help='map_collisions.<stratum>.tsv')
    p.add_argument('--map', required=True,
                   help='map.<stratum>.tsv.gz; read ONLY for REF lengths, to predict '
                        'rvtest\'s span-overlap over-count exactly')
    p.add_argument('--cohort', required=True)
    p.add_argument('--stratum', required=True)
    p.add_argument('--strict', dest='strict', action='store_true', default=True,
                   help='abort on any unexplained or missing identity row (default)')
    p.add_argument('--no-strict', dest='strict', action='store_false',
                   help='write the tables and report the failures without aborting')
    p.add_argument('--out-scan', default='gene_scan.tsv')
    p.add_argument('--out-qc', default='scan_qc.tsv')
    p.add_argument('--out-identity', default='map_identity.tsv')
    return p.parse_args()


def chrom_key(c):
    s = str(c)
    if s.isdigit():
        return (0, int(s), '')
    if s in ('', NA):
        return (2, 0, '')
    return (1, 0, s)


def warn(msg):
    sys.stdout.flush()
    print(msg, file=sys.stderr, flush=True)


def sample(items, limit=MAX_REPORT):
    s = sorted(str(x) for x in items)
    more = f' ... (+{len(s) - limit} more)' if len(s) > limit else ''
    return ', '.join(s[:limit]) + more


def parse_int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def parse_pvalue(v):
    """A finite p in [0, 1], or None. 'Returned a p' means this, not 'has a row':
    rvtest writes NA for a set it could not evaluate and still emits the row."""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    if x != x or x in (float('inf'), float('-inf')) or not 0.0 <= x <= 1.0:
        return None
    return x


def chi2_1df(p):
    """1-df chi-square quantile of a two-sided p, through the normal."""
    return NormalDist().inv_cdf(max(p, P_FLOOR) / 2.0) ** 2


def lambda_gc(ps):
    return statistics.median(chi2_1df(p) for p in ps) / CHI2_MEDIAN_1DF if ps else None


def split_fields(path, line, ln, n_expected):
    f = line.rstrip('\n').split('\t')
    if len(f) != n_expected:
        raise SystemExit(f'ABORT: {path}:{ln} has {len(f)} field(s), expected {n_expected}')
    return f


def read_gene_index(path):
    """set_name -> {gene_symbol, chrom, n_var_map, pos_min, pos_max}. This is the
    designed experiment: every count in scan_qc.tsv is relative to it."""
    idx = {}
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        col = {n: i for i, n in enumerate(head)}
        missing = [c for c in INDEX_COLUMNS if c not in col]
        if missing:
            raise SystemExit(f'ABORT: {path} has no {missing} column(s); found {head}')
        for ln, raw in enumerate(fh, 2):
            if not raw.strip():
                continue
            f = split_fields(path, raw, ln, len(head))
            sn = f[col['set_name']]
            if sn in idx:
                raise SystemExit(f'ABORT: {path}:{ln} lists set {sn!r} twice; the index is the '
                                 f'denominator, so a set counted twice inflates it')
            n_var = parse_int(f[col['n_var_map']])
            if n_var is None:
                raise SystemExit(f'ABORT: {path}:{ln} set {sn!r} has n_var_map '
                                 f'{f[col["n_var_map"]]!r}, which is not an integer')
            idx[sn] = {'gene_symbol': f[col['gene_symbol']], 'chrom': f[col['chrom']],
                       'n_var_map': n_var, 'pos_min': f[col['pos_min']],
                       'pos_max': f[col['pos_max']]}
    if not idx:
        raise SystemExit(f'ABORT: {path} lists no set; there is no experiment to scan')
    return idx


def read_denominator(path, cohort, stratum, n_index):
    """The one data row of denominator.<stratum>.tsv, checked against the index."""
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        if tuple(head) != DENOM_COLUMNS:
            raise SystemExit(f'ABORT: {path} header is {head}, expected {list(DENOM_COLUMNS)}')
        rows = [split_fields(path, raw, ln, len(head))
                for ln, raw in enumerate(fh, 2) if raw.strip()]
    if len(rows) != 1:
        raise SystemExit(f'ABORT: {path} has {len(rows)} data row(s); the denominator file '
                         f'carries exactly one, for exactly one (cohort, stratum)')
    row = dict(zip(head, rows[0]))
    if (row['cohort'], row['stratum']) != (cohort, stratum):
        raise SystemExit(f'ABORT: {path} is {row["cohort"]}/{row["stratum"]} but this task is '
                         f'{cohort}/{stratum}. The threshold would be the wrong stratum\'s.')
    n_mapped = parse_int(row['n_genes_mapped'])
    if n_mapped is None:
        raise SystemExit(f'ABORT: {path} n_genes_mapped is {row["n_genes_mapped"]!r}')
    if n_mapped != n_index:
        raise SystemExit(
            f'ABORT: {path} says n_genes_mapped={n_mapped:,} but the gene index lists '
            f'{n_index:,} set(s). Both are written by one MERGE_MAP task from one filtered '
            f'map, so these files are not from the same run and alpha/n_genes_mapped does '
            f'not correct the experiment that was actually tested.')
    return row


def read_collisions(path, cohort, stratum):
    """set_name -> the largest NumVar inflation a position shared by two GENES can
    produce: at most (n_variant_ids - 1) foreign records per colliding position."""
    bound = {}
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        if tuple(head) != COLLISION_COLUMNS:
            raise SystemExit(f'ABORT: {path} header is {head}, '
                             f'expected {list(COLLISION_COLUMNS)}')
        col = {n: i for i, n in enumerate(head)}
        n_pos = 0
        for ln, raw in enumerate(fh, 2):
            if not raw.strip():
                continue
            f = split_fields(path, raw, ln, len(head))
            names = [s for s in f[col['set_names']].split(',') if s]
            vids = [v for v in f[col['variant_ids']].split(',') if v]
            if len(names) < 2:
                raise SystemExit(f'ABORT: {path}:{ln} lists {len(names)} set name(s); this '
                                 f'file holds only positions carrying more than one')
            n_pos += 1
            for sn in names:
                bound[sn] = bound.get(sn, 0) + max(len(vids) - 1, 0)
    print(f'[scan] {cohort}/{stratum}: {n_pos:,} colliding position(s) touching '
          f'{len(bound):,} set(s) -> at most {sum(bound.values()):,} foreign record(s) '
          f'reachable by an rvtest range')
    return bound


def span_overlap_excess(path):
    """set_name -> the EXACT number of extra records rvtest will report in NumVar.

    excess = the number of records whose REF is longer than one base and whose
    span [p, p+L-1] STRICTLY contains another mapped variant's position. Strictness
    matters: two records at the SAME position (a multi-allelic site) do not
    over-count, because rvtest consolidates identical ranges. Validated against
    every rvtest result on the reference data: 81,746 / 81,746 genes exact, zero
    counter-examples.
    """
    opener = gzip.open if str(path).endswith('.gz') else open
    per_set = {}
    with opener(path, 'rt') as fh:
        head = fh.readline().rstrip('\n').split('\t')
        col = {n: i for i, n in enumerate(head)}
        for c in ('set_name', 'pos', 'variant_id'):
            if c not in col:
                raise SystemExit(f'ABORT: {path} has no {c} column; found {head}')
        i_sn, i_pos, i_vid = col['set_name'], col['pos'], col['variant_id']
        for ln, raw in enumerate(fh, 2):
            line = raw.rstrip('\n')
            if not line:
                continue
            f = line.split('\t')
            parts = f[i_vid].split(':')
            if len(parts) < 4:
                raise SystemExit(f'ABORT: {path}:{ln} variant_id {f[i_vid]!r} is not '
                                 f'chr<N>:<pos>:<REF>:<ALT>')
            try:
                pos = int(f[i_pos])
            except ValueError:
                raise SystemExit(f'ABORT: {path}:{ln} pos is {f[i_pos]!r}, not an integer')
            per_set.setdefault(f[i_sn], []).append((pos, len(parts[2])))
    out = {}
    for sn, rows in per_set.items():
        positions = {p for p, _L in rows}
        n = sum(1 for p, L in rows
                if L > 1 and any(p < q <= p + L - 1 for q in positions))
        if n:
            out[sn] = n
    return out


def read_scan(path, cohort, stratum):
    """One merged.tsv, on the gene_scan schema."""
    rows = []
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        if tuple(head) != SCAN_COLUMNS:
            raise SystemExit(
                f'ABORT: {path} header is {head}, expected the gene_scan schema '
                f'{list(SCAN_COLUMNS)}. This is a merged.tsv from MERGE_RVTEST, not a raw '
                f'engine output.')
        for ln, raw in enumerate(fh, 2):
            if not raw.strip():
                continue
            rows.append(dict(zip(head, split_fields(path, raw, ln, len(head)))))
    if not rows:
        raise SystemExit(f'ABORT: {path} carries no data row.')
    methods = {r['method'] for r in rows}
    if len(methods) != 1:
        raise SystemExit(f'ABORT: {path} carries methods {sample(methods)}; one merged.tsv is '
                         f'one method.')
    bad = methods - set(METHODS)
    if bad:
        raise SystemExit(f'ABORT: {path} method {sample(bad)} is not one of {list(METHODS)}')
    wrong = {(r['cohort'], r['stratum']) for r in rows} - {(cohort, stratum)}
    if wrong:
        raise SystemExit(f'ABORT: {path} carries {sample(wrong)} but this task is '
                         f'{cohort}/{stratum}')
    return rows


def classify(gi, n_map, n_rv, has_rv, bound, excess):
    if gi is None:
        # rvtest cannot invent a set name, so this came from a set file that a
        # different EMIT_TEST_INPUTS task wrote. There is no n_var_map to test
        # against and the denominator never counted it.
        return 'unexplained'
    if not has_rv:
        return 'missing'
    if n_rv is not None and n_rv == n_map:
        return 'ok'
    # The span-overlap predictor is an EQUALITY, so it is tested before the
    # collision bound, which is only an upper bound.
    if excess and n_rv is not None and n_rv == n_map + excess:
        return 'explained_by_span_overlap'
    if n_rv is not None and n_rv > n_map and bound is not None and n_rv - n_map <= bound:
        return 'explained_by_collision'
    return 'unexplained'


def build_identity(index, rows, bounds, excess, cohort, stratum):
    """One row per set in the union of the index and rvtest's output. The two
    methods read ONE set file, so their NumVar must agree; a disagreement is
    displayed as a/b and filed as unexplained rather than collapsed to a number."""
    counts, chroms = {}, {}
    for r in rows:
        counts.setdefault(r['set_name'], set()).add(r['n_var_engine'])
        chroms.setdefault(r['set_name'], set()).add(r['chrom'])
    universe = set(index) | set(counts)
    out, conflicts = [], []
    for sn in universe:
        gi = index.get(sn)
        n_map = gi['n_var_map'] if gi else None
        vals = sorted(counts.get(sn, set()))
        has_rv = sn in counts
        if len(vals) == 1:
            disp, n_rv = vals[0], parse_int(vals[0])
        elif vals:
            disp, n_rv = '/'.join(vals), None
            conflicts.append(f'{sn}={disp}')
        else:
            disp, n_rv = NA, None
        ex = excess.get(sn, 0)
        reason = classify(gi, n_map, n_rv, has_rv, bounds.get(sn), ex)
        chrom = gi['chrom'] if gi else (chroms[sn].pop() if len(chroms.get(sn, ())) == 1 else NA)
        out.append({
            'cohort': cohort, 'stratum': stratum, 'chrom': chrom, 'set_name': sn,
            'n_var_map': str(n_map) if n_map is not None else NA,
            'n_var_engine': disp, 'span_overlap_excess': str(ex),
            'ok': '1' if reason == 'ok' else '0', 'reason': reason,
            '_pos': parse_int(gi['pos_min']) if gi else None,
        })
    out.sort(key=lambda r: (chrom_key(r['chrom']),
                            r['_pos'] if r['_pos'] is not None else -1, r['set_name']))
    return out, conflicts


def bh_adjust(items, n):
    """Benjamini-Hochberg adjusted p-values (q-values), monotone, pure Python.

    n IS n_genes_mapped -- the family Bonferroni corrects over -- and the
    n - len(items) genes rvtest did not return enter the family at p = 1.0. They
    are not materialised, and they need not be: every imputed value is 1.0, which
    is >= any real p, so a returned gene at rank i among the returned keeps global
    rank i; for a tail position j > len(items), ranked_j = 1.0 * N / j >= 1.0, and
    every ranked value is clipped at 1.0, so the reverse running minimum arriving
    from the tail is exactly 1.0 and can never lower a returned gene's q.
    Computing q over the returned genes with multiplier N is therefore IDENTICAL
    to materialising the ones.

    items: {key: (p_float, p_str)}.  Returns {key: (q_float, q_str)}.
    """
    if not items or n <= 0:
        return {}
    order = sorted(items.items(), key=lambda kv: kv[1][0])   # stable
    ranked = [pv[0] * n / (i + 1) for i, (_k, pv) in enumerate(order)]
    run, out = 1.0, {}
    for i in range(len(ranked) - 1, -1, -1):
        run = min(run, ranked[i], 1.0)
        out[order[i][0]] = (run, f'{run:.12g}')
    return out


def build_qc(rows, index, denom, cohort, stratum):
    """One row per method. The population is the INDEXED sets: a set rvtest
    returned but the map never assigned cannot be counted as returned without
    letting it cancel a set that really is missing."""
    thr_str = denom['bonferroni_threshold']
    thr = None if thr_str == NA else float(thr_str)
    n_mapped = int(denom['n_genes_mapped'])
    alpha = float(denom['alpha'])
    if not 0.0 < alpha < 1.0:
        raise SystemExit(f'ABORT: denominator alpha is {alpha!r}; expected 0 < alpha < 1')
    if thr is not None and abs(alpha / n_mapped - thr) > REL_TOL * max(thr, 1e-300):
        raise SystemExit(f'ABORT: denominator says alpha {alpha} over {n_mapped} gene(s) but '
                         f'bonferroni_threshold {thr_str}; {alpha / n_mapped:.12g} expected.')

    groups = {}
    for r in rows:
        g = groups.setdefault(r['method'], {})
        if r['set_name'] not in index:
            continue
        p = parse_pvalue(r['pvalue'])
        if p is not None:
            g[r['set_name']] = (p, r['pvalue'])

    qc, bh = [], {}
    for method, ps in groups.items():
        n_ret = len(ps)
        if ps:
            lam = lambda_gc([p for p, _s in ps.values()])
            lam2 = lambda_gc([p for k, (p, _s) in ps.items() if index[k]['n_var_map'] == 2])
            lam3 = lambda_gc([p for k, (p, _s) in ps.items() if index[k]['n_var_map'] >= 3])
            # min_p is emitted VERBATIM as the engine wrote it.
            min_p = min(ps.values(), key=lambda v: v[0])[1]
            n_sig = sum(1 for p, _s in ps.values() if thr is not None and p < thr)
            qs = bh_adjust(ps, n_mapped)
            rejected = [ps[k][0] for k, (q, _qs) in qs.items() if q < alpha]
            n_sig_bh = len(rejected)
            # The LARGEST p BH rejected: the cut a Manhattan panel draws its
            # second line at; NA when nothing was rejected.
            if rejected:
                top = max(rejected)
                bh_cut = next(sv for pv, sv in ps.values() if pv == top)
            else:
                bh_cut = NA
        else:
            lam = lam2 = lam3 = None
            min_p, n_sig, qs, n_sig_bh, bh_cut = NA, 0, {}, 0, NA
        fmt = lambda x: f'{x:.4f}' if x is not None else NA
        qc.append({
            'cohort': cohort, 'stratum': stratum, 'method': method,
            'n_genes_mapped': str(n_mapped),
            'n_genes_returned': str(n_ret), 'n_genes_missing': str(n_mapped - n_ret),
            'lambda_gc': fmt(lam), 'lambda_gc_nvar2': fmt(lam2), 'lambda_gc_nvar_ge3': fmt(lam3),
            'min_p': min_p, 'n_significant_bonferroni': str(n_sig),
            'n_significant_bh': str(n_sig_bh), 'bh_threshold': bh_cut,
            'bonferroni_threshold': thr_str, 'min_num_var': denom['min_num_var'],
        })
        bh[method] = qs
    qc.sort(key=lambda r: METHODS.index(r['method']))
    return qc, thr, alpha, bh


def write_table(path, columns, rows):
    with open(path, 'w') as fh:
        fh.write('\t'.join(columns) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in columns) + '\n')


def main():
    args = parse_args()
    index = read_gene_index(args.gene_index)
    denom = read_denominator(args.denominator, args.cohort, args.stratum, len(index))
    bounds = read_collisions(args.collisions, args.cohort, args.stratum)

    rows, origin = [], {}
    for path in args.assoc:
        part = read_scan(path, args.cohort, args.stratum)
        for r in part:
            key = (r['set_name'], r['method'])
            if key in origin:
                raise SystemExit(f'ABORT: {key[0]} / {key[1]} is returned by both {origin[key]} '
                                 f'and {path}. gene_scan.tsv is one row per gene x method.')
            origin[key] = path
        rows.extend(part)
        print(f'[scan] {path}: {part[0]["method"]}, {len(part):,} row(s)')
    methods_seen = {r['method'] for r in rows}
    if methods_seen != set(METHODS):
        raise SystemExit(f'ABORT: expected one merged.tsv per method {list(METHODS)}, got '
                         f'{sorted(methods_seen)}')

    # ---- the map identity, written before anything else can abort on it.
    excess = span_overlap_excess(args.map)
    if excess:
        print(f'[scan] {len(excess):,} set(s) carry a multi-base-REF record whose span covers '
              f'another mapped variant; rvtest over-counts those by '
              f'{sum(excess.values()):,} record(s) in total (exact, see span_overlap_excess)')
    identity, conflicts = build_identity(index, rows, bounds, excess, args.cohort, args.stratum)
    write_table(args.out_identity, IDENTITY_COLUMNS, identity)
    tally = {k: 0 for k in REASONS}
    for r in identity:
        tally[r['reason']] += 1
    print(f'[scan] map identity over {len(identity):,} set(s): '
          + ' '.join(f'{k}={tally[k]:,}' for k in REASONS) + f' -> {args.out_identity}')
    if conflicts:
        warn(f'[scan] WARNING: {len(conflicts):,} set(s) get different variant counts from '
             f'the two methods, which read one set file: {sample(conflicts)}')

    failures = [r for r in identity if r['reason'] in FAIL_REASONS]
    if failures and args.strict:
        lines = [f'  chr{r["chrom"]} {r["set_name"]}: n_var_map={r["n_var_map"]} '
                 f'n_var_engine={r["n_var_engine"]} [{r["reason"]}]'
                 for r in failures[:MAX_REPORT]]
        more = (f'\n  ... (+{len(failures) - MAX_REPORT} more; every set is in '
                f'{args.out_identity})' if len(failures) > MAX_REPORT else '')
        raise SystemExit(
            f'ABORT: the map identity does not hold for {len(failures):,} of '
            f'{len(identity):,} set(s) in {args.cohort}/{args.stratum}.\n'
            + '\n'.join(lines) + more + '\n'
            f'       {args.out_identity} is written with every row, passing or not. '
            f'{args.out_scan} is NOT: a scan whose tested sets are not proven to be the '
            f'mapped sets is not a result. Re-run with --no-strict to emit it anyway.')
    if failures:
        warn(f'[scan] WARNING: --no-strict, so {len(failures):,} identity failure(s) are '
             f'reported and the scan is written anyway')

    # ---- the scan and its QC.
    qc, thr, alpha, bh = build_qc(rows, index, denom, args.cohort, args.stratum)
    thr_str = denom['bonferroni_threshold']
    for r in rows:
        p = parse_pvalue(r['pvalue'])
        r['bonferroni_threshold'] = thr_str
        # NA, not 0, when there is nothing to compare: a set rvtest could not
        # evaluate is already counted in n_genes_missing, and a 0 here would read
        # as 'tested, did not clear the bar'.
        r['significant_bonferroni'] = (
            NA if (p is None or thr is None) else ('1' if p < thr else '0'))
        q = bh.get(r['method'], {}).get(r['set_name'])
        r['fdr_bh'] = NA if q is None else q[1]
        r['significant_bh'] = NA if q is None else ('1' if q[0] < alpha else '0')

    # Both rules share alpha and n, so BH's smallest critical value is alpha/n --
    # exactly the Bonferroni threshold -- and every Bonferroni call is a BH call.
    # The containment is an identity; a violation can only mean the two rules saw
    # different families or different alphas.
    bad = [r for r in rows if r['significant_bonferroni'] == '1' and r['significant_bh'] != '1']
    if bad:
        named = [f'{r["method"]} {r["set_name"]} p={r["pvalue"]} q={r["fdr_bh"]}' for r in bad]
        raise SystemExit(
            f'ABORT: {len(bad):,} row(s) are Bonferroni-significant but not BH-significant, '
            f'which is impossible when both rules use the same alpha and family size. '
            f'{sample(named)}')

    def scan_key(r):
        pos = parse_int(r['pos_min'])
        return (chrom_key(r['chrom']), pos if pos is not None else -1, r['set_name'],
                METHODS.index(r['method']))

    rows.sort(key=scan_key)
    write_table(args.out_scan, SCAN_COLUMNS, rows)
    write_table(args.out_qc, QC_COLUMNS, qc)

    print(f'[scan] {args.cohort}/{args.stratum}: {len(rows):,} row(s) over {len(qc)} method(s), '
          f'judged by Bonferroni p < {thr_str} (alpha {denom["alpha"]} / '
          f'{denom["n_genes_mapped"]} gene(s)) and BH q < {denom["alpha"]} over the same '
          f'family -> {args.out_scan}')
    for r in qc:
        print(f'    {r["method"]:<6} returned {r["n_genes_returned"]:>6} / '
              f'{r["n_genes_mapped"]:<6} missing {r["n_genes_missing"]:>4}  '
              f'lambda {r["lambda_gc"]:>7} (n_var=2 {r["lambda_gc_nvar2"]}, >=3 '
              f'{r["lambda_gc_nvar_ge3"]})  min_p {r["min_p"]:<12} '
              f'bonferroni {r["n_significant_bonferroni"]:>4}  bh {r["n_significant_bh"]:>4} '
              f'(cut {r["bh_threshold"]})')
    n_missing = sum(int(r['n_genes_missing']) for r in qc)
    if n_missing:
        warn(f'[scan] NOTE: {n_missing:,} (method, gene) reading(s) returned no usable p. '
             f'n_genes_mapped stays the denominator; the shortfall is in {args.out_qc}.')
    print(f'[scan] {len(qc)} QC row(s) -> {args.out_qc}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in scan: {e}', file=sys.stderr)
        sys.exit(1)
