#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The deliverable. Every gene x cohort x stratum x method in ONE
#           evidence table; a robustness TIER per (gene, stratum); and, for the
#           tiered genes, the publication summary table -- one row per gene, one
#           block per stratum, in the layout the manuscript uses.
#
#           WHAT "ROBUST" MEANS HERE, AND WHAT IT DOES NOT. The scan already paid
#           alpha / n_genes_mapped per cohort x stratum and called a gene when
#           BH q < alpha (Bonferroni is reported as the stricter tier, and never
#           enters the rule). This step introduces NO new threshold and NO new
#           p-value. It asks one further question: does the call survive a
#           change of STATISTIC (CMC burden vs SKAT-O) and a change of SAMPLE
#           SELECTION (the three cohorts)? The tiers, from tier_of():
#
#             Tier 1  called in EVERY cohort, and in at least one cohort by
#                     BOTH methods
#             Tier 2  called in EVERY cohort, by at least one method each
#             Tier 3  called in at least two cohorts
#
#           THE COHORTS ARE NESTED (narrow c intermediate c full, differing
#           almost entirely in controls), so "called in every cohort" means the
#           signal survives ADDING CONTROLS. That is a sample-selection
#           sensitivity analysis, and it is not replication; every table this
#           step writes says so in its own caption.
#
#           not_in_map IS A STATE, NOT AN ABSENCE. The gene sets are not nested
#           across cohorts: a gene can have >= MinNumVar variants in one cohort's
#           map and fewer in another's (the minAC floor is per cohort). A gene
#           never tested in a cell must not read as "tested and not significant"
#           there, so the state is decided from gene_index_all -- the designed
#           experiment -- and never from a missing scan row.
#
#           A SAMPLE CARRYING TWO SITES OF ONE GENE. CMC collapses a sample to
#           carrier / non-carrier, so such a sample counts ONCE in the carrier
#           column and TWICE in the cumulative MAC. Where that happens in a
#           tiered gene the summary table says so in a column of its own, and
#           multi_carrier_detail.tsv lists the samples and their variants.
#
#           EVERY NUMBER IN THE SUMMARY TABLE HAS ONE SOURCE.
#             n_variants          gene_index_all n_var_map
#             carriers, MAC       set_stats via gene_scan (EMIT_TEST_INPUTS bed
#                                 decode, cross-checked per variant against
#                                 plink2 --freq counts)
#             MAF                 MAC / (2 * N_group), N from cohort_samples.tsv
#                                 -- a FIXED denominator, not OBS_CT, so the two
#                                 groups' MAFs are on the same footing and the
#                                 ~2 % genotype missingness is not a per-gene
#                                 wobble in the third decimal
#             SKAT-O rho, p       gene_scan skato row (rvtest's own rho grid)
#             CMC OR (95 % CI)    exp(beta), exp(beta -/+ 1.96 se) from rvtest's
#                                 --burden cmcWald on the SAME collapsed genotype
#             CMC p               the CMC SCORE p (not the Wald p; the two differ)
# Component: assoc_gene
# ---------------------------------------------------------------------------
import argparse
import math
import sys
from collections import OrderedDict, defaultdict

METHODS = ('cmc', 'skato')
STATES = ('called_bonferroni', 'called_bh_only', 'not_significant', 'not_in_map')
CALLED_STATES = ('called_bonferroni', 'called_bh_only')

SCAN_COLUMNS = ('cohort', 'stratum', 'method', 'set_name', 'gene_symbol', 'chrom',
                'pos_min', 'pos_max', 'n_var_map', 'n_var_engine', 'n_informative',
                'pvalue', 'rho', 'beta', 'se', 'or', 'or_l95', 'or_u95',
                'mac', 'mac_case', 'mac_control', 'n_carrier_case', 'n_carrier_control',
                'bonferroni_threshold', 'significant_bonferroni', 'fdr_bh', 'significant_bh')
QC_READ = ('cohort', 'stratum', 'method', 'lambda_gc', 'n_genes_mapped')
INDEX_READ = ('cohort', 'stratum', 'set_name', 'gene_symbol', 'chrom', 'n_var_map', 'pos_min')
SAMPLES_COLUMNS = ('cohort', 'n_samples', 'n_case', 'n_control')
DEN_READ = ('cohort', 'stratum', 'n_genes_mapped', 'alpha', 'bonferroni_threshold')
SET_STATS_READ = ('cohort', 'stratum', 'set_name', 'n_multi_case', 'n_multi_control',
                  'max_var_per_sample_case', 'max_var_per_sample_control')
WINDOW_COLUMNS = ('gene_symbol', 'set_name', 'chrom', 'start', 'end')
# The gene model panel needs the gene's own extent, which the variant span only
# bounds from the inside; a pad of this size around it reaches the promoter and
# the last exon of every gene in the current tier list without pulling in a
# neighbourhood.
WINDOW_PAD_FRAC, WINDOW_PAD_MIN = 0.10, 20_000
VARIANT_KEY = ('cohort', 'stratum', 'set_name')
MULTI_KEY = ('cohort', 'stratum', 'set_name')

EVIDENCE_COLUMNS = ('gene_symbol', 'set_name', 'chrom', 'cohort', 'stratum', 'method', 'state',
                    'pvalue', 'nlp', 'fdr_bh', 'n_var_map', 'mac_case', 'mac_control',
                    'n_carrier_case', 'n_carrier_control', 'rho', 'or', 'or_l95', 'or_u95',
                    'lambda_gc')

TIER_COLUMNS = ('gene_symbol', 'set_name', 'chrom', 'stratum', 'tier', 'n_cohorts_called',
                'cohorts_called', 'n_cohorts_called_cmc', 'n_cohorts_called_skato',
                'n_cohorts_both_methods', 'called_bonferroni_anywhere', 'n_cohorts_mapped',
                'cohorts_not_in_map', 'min_p_cmc', 'min_p_skato', 'best_cohort',
                'n_multi_site_carriers_case_max', 'n_multi_site_carriers_control_max', 'note')

# One block per stratum in the summary table, in this order. The first twelve
# are the manuscript's columns; the last is the multi-site carrier count that
# explains any gap between carriers and MAC.
BLOCK_COLUMNS = ('tier', 'n_variants', 'carriers_case', 'carriers_control', 'mac_case',
                 'mac_control', 'maf_case', 'maf_control', 'skato_rho', 'skato_p',
                 'cmc_or_95ci', 'cmc_p', 'multi_site_carriers')
BLOCK_HEADERS = ('Tier', 'Variants', 'Carriers case (%)', 'Carriers control (%)', 'MAC case',
                 'MAC control', 'MAF case', 'MAF control', 'SKAT-O rho', 'SKAT-O P',
                 'CMC OR (95% CI)', 'CMC P', 'Multi-site carriers case / control')
SUMMARY_TIERS = (1, 2)

NA = 'NA'
DASH = '—'
P_FLOOR = 1e-300
REL_TOL = 1e-9


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gene-scan-all', required=True)
    p.add_argument('--scan-qc-all', required=True)
    p.add_argument('--gene-index-all', required=True,
                   help='gene_index_all.tsv, keyed by cohort and stratum; the designed experiment')
    p.add_argument('--cohort-samples', required=True,
                   help='cohort_samples.tsv from COLLECT_ALL; N_case / N_control per cohort')
    p.add_argument('--denominators', required=True,
                   help='denominators.tsv; alpha is asserted against --alpha')
    p.add_argument('--set-stats-all', default=None,
                   help='set_stats_all.tsv; multi-site carrier counts per set (optional only '
                        'so an older table set still runs; the summary column is then NA)')
    p.add_argument('--variant-stats-all', required=True)
    p.add_argument('--multi-carriers-all', required=True)
    p.add_argument('--cohorts', nargs='+', required=True,
                   help='cohort tags in the order the summary blocks and tier rule use')
    p.add_argument('--strata', nargs='+', required=True, help='stratum tags, widest first')
    p.add_argument('--report-cohort', required=True,
                   help='the cohort whose numbers fill summary_table.md; every cohort gets '
                        'its own .tsv regardless')
    p.add_argument('--alpha', type=float, default=0.05)
    p.add_argument('--out-evidence', default='evidence.tsv')
    p.add_argument('--out-tiers', default='tiers.tsv')
    p.add_argument('--out-variant-detail', default='variant_detail.tsv')
    p.add_argument('--out-multi-detail', default='multi_carrier_detail.tsv')
    p.add_argument('--out-gene-windows', default='gene_windows.tsv',
                   help='one row per Tier 1/2 gene: the window whose Ensembl gene models the '
                        'per-gene figure draws under its variant panel')
    p.add_argument('--out-summary-prefix', default='summary_table')
    return p.parse_args()


# ── the rule ─────────────────────────────────────────────────────────────────
def tier_of(called, n_cohorts_total):
    """The robustness tier of one (gene, stratum).

    called: {cohort: {'cmc': bool, 'skato': bool}} over the cohorts in which the
    gene was MAPPED (a cohort where it is not in the map is absent from the
    dict, and therefore cannot count as called). n_cohorts_total is the number
    of cohorts in the run, so Tier 1 and Tier 2 require a call in every one.

    Pure, so verify.sh can fixture-test it without a scan on disk.
    """
    any_c = [c for c, m in called.items() if m['cmc'] or m['skato']]
    both_c = [c for c, m in called.items() if m['cmc'] and m['skato']]
    if len(any_c) == n_cohorts_total and both_c:
        return 1
    if len(any_c) == n_cohorts_total:
        return 2
    if len(any_c) >= 2:
        return 3
    return 0


# ── readers ──────────────────────────────────────────────────────────────────
def read_tsv(path, required, label):
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        missing = [c for c in required if c not in head]
        if missing:
            raise SystemExit(f'ABORT: {path} has no {missing} column(s); found {head}. '
                             f'This is not {label}.')
        rows = []
        for ln, raw in enumerate(fh, 2):
            line = raw.rstrip('\n')
            if not line:
                continue
            f = line.split('\t')
            if len(f) != len(head):
                raise SystemExit(f'ABORT: {path}:{ln} has {len(f)} field(s), expected '
                                 f'{len(head)}')
            rows.append(dict(zip(head, f)))
    return rows, head


def as_float(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if x != x or x in (float('inf'), float('-inf')) else x


def as_int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def chrom_key(c):
    s = str(c)
    if s.isdigit():
        return (0, int(s), '')
    if s in ('', NA):
        return (2, 0, '')
    return (1, 0, s)


def sample(items, limit=12):
    s = sorted(str(x) for x in items)
    more = f' ... (+{len(s) - limit} more)' if len(s) > limit else ''
    return ', '.join(s[:limit]) + more


def write_table(path, columns, rows):
    with open(path, 'w') as fh:
        fh.write('\t'.join(columns) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in columns) + '\n')


# ── formats: the manuscript's ────────────────────────────────────────────────
def fmt_carriers(n, n_group):
    n = as_int(n)
    if n is None or not n_group:
        return DASH
    return f'{n} ({100.0 * n / n_group:.2f}%)'


def fmt_maf(mac, n_group):
    m = as_int(mac)
    if m is None or not n_group:
        return DASH
    return f'{m / (2.0 * n_group):.3f}'


def fmt_p(p):
    x = as_float(p)
    return DASH if x is None else f'{x:.2E}'


def fmt_or(o, lo, hi):
    o, lo, hi = as_float(o), as_float(lo), as_float(hi)
    if o is None or lo is None or hi is None:
        return DASH
    return f'{o:.2f} ({lo:.2f}-{hi:.2f})'


def fmt_int(v):
    x = as_int(v)
    return DASH if x is None else str(x)


def fmt_rho(v):
    x = as_float(v)
    return DASH if x is None else f'{x:g}'


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    cohorts, strata = list(args.cohorts), list(args.strata)
    if len(set(cohorts)) != len(cohorts) or len(set(strata)) != len(strata):
        raise SystemExit('ABORT: --cohorts / --strata repeat a tag')
    if args.report_cohort not in cohorts:
        raise SystemExit(f'ABORT: --report-cohort {args.report_cohort!r} is not one of '
                         f'{cohorts}')
    if not 0.0 < args.alpha < 1.0:
        raise SystemExit(f'ABORT: --alpha must lie in (0, 1); got {args.alpha}')

    # ---- the designed experiment: which gene is mapped in which cell.
    idx_rows, _h = read_tsv(args.gene_index_all, INDEX_READ, 'gene_index_all.tsv')
    mapped = defaultdict(dict)          # (cohort, stratum) -> set_name -> row
    gene_info = {}                      # (stratum, set_name) -> (symbol, chrom, pos_min)
    for r in idx_rows:
        key = (r['cohort'], r['stratum'])
        if r['cohort'] not in cohorts or r['stratum'] not in strata:
            continue
        if r['set_name'] in mapped[key]:
            raise SystemExit(f'ABORT: {args.gene_index_all} lists {r["set_name"]!r} twice in '
                             f'{key[0]}/{key[1]}')
        mapped[key][r['set_name']] = r
        gk = (r['stratum'], r['set_name'])
        info = (r['gene_symbol'], r['chrom'], as_int(r['pos_min']) or 0)
        if gk in gene_info and gene_info[gk][:2] != info[:2]:
            raise SystemExit(f'ABORT: set {r["set_name"]!r} is {gene_info[gk][:2]} in one '
                             f'cohort and {info[:2]} in another; a set name must mean one '
                             f'gene on one chromosome everywhere')
        gene_info.setdefault(gk, info)
    for c in cohorts:
        for s in strata:
            if not mapped.get((c, s)):
                raise SystemExit(f'ABORT: {args.gene_index_all} has no gene for {c}/{s}; the '
                                 f'evidence table cannot be complete')
    universe = {s: sorted({sn for (c, st), d in mapped.items() if st == s for sn in d},
                          key=lambda sn: (chrom_key(gene_info[(s, sn)][1]),
                                          gene_info[(s, sn)][2], sn))
                for s in strata}

    # ---- sample counts and alpha.
    smp_rows, _h = read_tsv(args.cohort_samples, SAMPLES_COLUMNS, 'cohort_samples.tsv')
    n_of = {r['cohort']: (as_int(r['n_case']), as_int(r['n_control']), as_int(r['n_samples']))
            for r in smp_rows}
    for c in cohorts:
        if c not in n_of or None in n_of[c] or n_of[c][0] + n_of[c][1] != n_of[c][2]:
            raise SystemExit(f'ABORT: {args.cohort_samples} has no usable row for {c!r}')
    den_rows, _h = read_tsv(args.denominators, DEN_READ, 'denominators.tsv')
    for r in den_rows:
        a = as_float(r['alpha'])
        if a is None or abs(a - args.alpha) > REL_TOL * args.alpha:
            raise SystemExit(f'ABORT: {args.denominators} {r["cohort"]}/{r["stratum"]} fixed '
                             f'alpha={r["alpha"]!r} but --alpha is {args.alpha:g}; the tier '
                             f'rule would describe calls made at another alpha')

    # ---- the scan and its calibration.
    qc_rows, _h = read_tsv(args.scan_qc_all, QC_READ, 'scan_qc_all.tsv')
    lam = {(r['cohort'], r['stratum'], r['method']): r['lambda_gc'] for r in qc_rows}
    scan_rows, head = read_tsv(args.gene_scan_all, SCAN_COLUMNS, 'gene_scan_all.tsv')
    if tuple(head) != SCAN_COLUMNS:
        raise SystemExit(f'ABORT: {args.gene_scan_all} header is {head}, expected the '
                         f'gene_scan schema {list(SCAN_COLUMNS)}')
    scan = {}
    for r in scan_rows:
        k = (r['cohort'], r['stratum'], r['method'], r['set_name'])
        if k in scan:
            raise SystemExit(f'ABORT: {args.gene_scan_all} carries {k} twice')
        scan[k] = r
        if r['method'] not in METHODS:
            raise SystemExit(f'ABORT: method {r["method"]!r} is not one of {METHODS}')
        cell = mapped.get((r['cohort'], r['stratum']))
        if cell is not None and r['set_name'] not in cell:
            raise SystemExit(f'ABORT: {r["cohort"]}/{r["stratum"]} scan returns '
                             f'{r["set_name"]!r}, which its gene index never mapped')

    # ---- per-set multi-site carrier counts (from set_stats, when given).
    multi_counts = {}
    if args.set_stats_all:
        ss_rows, _h = read_tsv(args.set_stats_all, SET_STATS_READ, 'set_stats_all.tsv')
        for r in ss_rows:
            multi_counts[(r['cohort'], r['stratum'], r['set_name'])] = (
                as_int(r['n_multi_case']), as_int(r['n_multi_control']),
                as_int(r['max_var_per_sample_case']), as_int(r['max_var_per_sample_control']))

    # ---- the evidence table: the full design, one row per cell x method.
    evidence, called = [], defaultdict(dict)     # (stratum, set_name) -> cohort -> {m: bool}
    bonf_any, min_p = defaultdict(bool), {}
    n_state = defaultdict(int)
    for s in strata:
        for sn in universe[s]:
            symbol, chrom, _pos = gene_info[(s, sn)]
            for c in cohorts:
                in_map = sn in mapped[(c, s)]
                if in_map:
                    called[(s, sn)][c] = {}
                for m in METHODS:
                    r = scan.get((c, s, m, sn))
                    if not in_map:
                        state = 'not_in_map'
                    elif r is None:
                        # Mapped, never returned: SCAN entered it at p = 1 in the
                        # BH family, so it is a tested-and-not-called gene.
                        state = 'not_significant'
                    elif r['significant_bonferroni'] == '1':
                        state = 'called_bonferroni'
                    elif r['significant_bh'] == '1':
                        state = 'called_bh_only'
                    else:
                        state = 'not_significant'
                    n_state[state] += 1
                    if in_map:
                        called[(s, sn)][c][m] = state in CALLED_STATES
                        if state == 'called_bonferroni':
                            bonf_any[(s, sn)] = True
                    p = as_float(r['pvalue']) if r else None
                    if p is not None:
                        k = (s, sn, m)
                        if k not in min_p or p < min_p[k][0]:
                            min_p[k] = (p, r['pvalue'], c)
                    g = lambda col: (r[col] if r else NA)
                    evidence.append({
                        'gene_symbol': symbol, 'set_name': sn, 'chrom': chrom,
                        'cohort': c, 'stratum': s, 'method': m, 'state': state,
                        'pvalue': g('pvalue'),
                        'nlp': f'{-math.log10(max(p, P_FLOOR)):.4f}' if p is not None else NA,
                        'fdr_bh': g('fdr_bh'),
                        'n_var_map': mapped[(c, s)][sn]['n_var_map'] if in_map else NA,
                        'mac_case': g('mac_case'), 'mac_control': g('mac_control'),
                        'n_carrier_case': g('n_carrier_case'),
                        'n_carrier_control': g('n_carrier_control'),
                        'rho': g('rho'), 'or': g('or'), 'or_l95': g('or_l95'),
                        'or_u95': g('or_u95'),
                        'lambda_gc': lam.get((c, s, m), NA) if in_map else NA,
                    })
    n_expected = sum(len(universe[s]) for s in strata) * len(cohorts) * len(METHODS)
    if len(evidence) != n_expected:
        raise SystemExit(f'ABORT: {len(evidence):,} evidence row(s) for a design of '
                         f'{n_expected:,}')
    write_table(args.out_evidence, EVIDENCE_COLUMNS, evidence)

    # ---- tiers.
    tiers, tier_of_gene = [], {}
    for s in strata:
        for sn in universe[s]:
            symbol, chrom, _pos = gene_info[(s, sn)]
            cl = called[(s, sn)]
            t = tier_of(cl, len(cohorts))
            tier_of_gene[(s, sn)] = t
            if t == 0:
                continue
            any_c = [c for c in cohorts if c in cl and (cl[c]['cmc'] or cl[c]['skato'])]
            both_c = [c for c in cohorts if c in cl and cl[c]['cmc'] and cl[c]['skato']]
            not_in = [c for c in cohorts if c not in cl]
            pc, ps = min_p.get((s, sn, 'cmc')), min_p.get((s, sn, 'skato'))
            best = min((x for x in (pc, ps) if x), key=lambda x: x[0], default=None)
            mc = [multi_counts.get((c, s, sn)) for c in cohorts]
            mc = [x for x in mc if x and None not in x]
            notes = []
            no_or = [c for c in cohorts if c in cl
                     and (scan.get((c, s, 'cmc', sn)) or {}).get('or', NA) == NA]
            if no_or:
                notes.append('cmc_or_na:' + '|'.join(no_or))
            if mc and max(x[0] + x[1] for x in mc):
                notes.append('multi_site_carriers')
            tiers.append({
                'gene_symbol': symbol, 'set_name': sn, 'chrom': chrom, 'stratum': s,
                'tier': t, 'n_cohorts_called': len(any_c), 'cohorts_called': ','.join(any_c),
                'n_cohorts_called_cmc': sum(1 for c in cl if cl[c]['cmc']),
                'n_cohorts_called_skato': sum(1 for c in cl if cl[c]['skato']),
                'n_cohorts_both_methods': len(both_c),
                'called_bonferroni_anywhere': int(bonf_any[(s, sn)]),
                'n_cohorts_mapped': len(cl), 'cohorts_not_in_map': ','.join(not_in),
                'min_p_cmc': pc[1] if pc else NA, 'min_p_skato': ps[1] if ps else NA,
                'best_cohort': best[2] if best else NA,
                'n_multi_site_carriers_case_max': max((x[0] for x in mc), default=NA),
                'n_multi_site_carriers_control_max': max((x[1] for x in mc), default=NA),
                'note': ';'.join(notes),
            })
    tiers.sort(key=lambda r: (r['tier'], strata.index(r['stratum']), chrom_key(r['chrom']),
                              gene_info[(r['stratum'], r['set_name'])][2], r['set_name']))
    write_table(args.out_tiers, TIER_COLUMNS, tiers)
    tiered = {(r['stratum'], r['set_name']) for r in tiers}

    # ---- the summary tables: Tier 1 + 2 genes, every cohort, ReportCohort as .md.
    summary_sets = OrderedDict()
    for r in tiers:
        if r['tier'] in SUMMARY_TIERS:
            summary_sets.setdefault(r['set_name'], r['gene_symbol'])
    # A gene is one row whatever stratum tiered it, so the row's chromosome and
    # position are taken from whichever stratum knows it.
    def gene_key(sn):
        for s in strata:
            if (s, sn) in gene_info:
                return (chrom_key(gene_info[(s, sn)][1]), gene_info[(s, sn)][2], sn)
        return ((3, 0, ''), 0, sn)
    best_tier = {sn: min(tier_of_gene.get((s, sn), 0) or 9 for s in strata)
                 for sn in summary_sets}
    order = sorted(summary_sets, key=lambda sn: (best_tier[sn], gene_key(sn)))

    def block(cohort, s, sn):
        n_case, n_ctrl, _n = n_of[cohort]
        if sn not in mapped[(cohort, s)]:
            return {c: DASH for c in BLOCK_COLUMNS}
        t = tier_of_gene.get((s, sn), 0)
        cmc = scan.get((cohort, s, 'cmc', sn)) or {}
        sk = scan.get((cohort, s, 'skato', sn)) or {}
        src = cmc or sk
        mc = multi_counts.get((cohort, s, sn))
        return {
            'tier': str(t),
            'n_variants': mapped[(cohort, s)][sn]['n_var_map'],
            'carriers_case': fmt_carriers(src.get('n_carrier_case'), n_case),
            'carriers_control': fmt_carriers(src.get('n_carrier_control'), n_ctrl),
            'mac_case': fmt_int(src.get('mac_case')),
            'mac_control': fmt_int(src.get('mac_control')),
            'maf_case': fmt_maf(src.get('mac_case'), n_case),
            'maf_control': fmt_maf(src.get('mac_control'), n_ctrl),
            'skato_rho': fmt_rho(sk.get('rho')),
            'skato_p': fmt_p(sk.get('pvalue')),
            'cmc_or_95ci': fmt_or(cmc.get('or'), cmc.get('or_l95'), cmc.get('or_u95')),
            'cmc_p': fmt_p(cmc.get('pvalue')),
            'multi_site_carriers': (f'{mc[0]} / {mc[1]}' if mc and None not in mc[:2]
                                    else NA),
        }

    columns = ['gene_symbol', 'chrom'] + [f'{s}.{c}' for s in strata for c in BLOCK_COLUMNS]
    tables = {}
    for cohort in cohorts:
        rows = []
        for sn in order:
            row = {'gene_symbol': summary_sets[sn],
                   'chrom': next(gene_info[(s, sn)][1] for s in strata if (s, sn) in gene_info)}
            for s in strata:
                for c, v in block(cohort, s, sn).items():
                    row[f'{s}.{c}'] = v
            rows.append(row)
        tables[cohort] = rows
        write_table(f'{args.out_summary_prefix}.{cohort}.tsv', columns, rows)

    # summary_table.md -- the ReportCohort, one markdown table per stratum so
    # the row does not run to 28 columns.
    rc = args.report_cohort
    n_case, n_ctrl, n_all = n_of[rc]
    with open(f'{args.out_summary_prefix}.md', 'w') as fh:
        fh.write(f'# Robust genes — {rc}\n\n')
        fh.write(f'{n_all:,} samples ({n_case:,} cases, {n_ctrl:,} controls). '
                 f'{len(order)} gene(s) at Tier 1 or Tier 2 in at least one stratum.\n\n')
        if not order:
            fh.write('_No gene reached Tier 1 or Tier 2._\n\n')
        for s in strata:
            fh.write(f'## {s.replace("_", " + ").upper()}\n\n')
            fh.write('| Gene | Chr | ' + ' | '.join(BLOCK_HEADERS) + ' |\n')
            fh.write('|' + '---|' * (2 + len(BLOCK_HEADERS)) + '\n')
            for row in tables[rc]:
                cells = [f'*{row["gene_symbol"]}*', row['chrom']]
                cells += [row[f'{s}.{c}'] for c in BLOCK_COLUMNS]
                fh.write('| ' + ' | '.join(cells) + ' |\n')
            fh.write('\n')
        fh.write('**Definitions.** Variants = variants of the stratum mapped to the gene in '
                 'this cohort (>= MinNumVar). Carriers = samples with >= 1 minor allele at '
                 'any of them (%: of the group). MAC = cumulative minor-allele count over '
                 'the set. MAF = MAC / (2 x N of the group), a fixed denominator. SKAT-O '
                 'rho = the mixing parameter rvtest selected; P = SKAT-O p. CMC OR (95% CI) '
                 '= exp(beta) from rvtest --burden cmcWald with a Wald normal CI; CMC P = '
                 'the CMC score-test p (the Wald p is not shown). Multi-site carriers = '
                 'samples carrying a minor allele at >= 2 distinct sites of the gene '
                 '(case / control); CMC counts each once, MAC counts every site. '
                 f'{DASH} = the gene is not in this cohort\'s map for that stratum, or the '
                 'value is undefined.\n\n')
        fh.write('**Tiers.** A gene is *called* in a cohort x stratum x method when '
                 f'Benjamini-Hochberg q < {args.alpha:g} over that cell\'s mapped genes. '
                 'Tier 1 = called in every cohort and by both methods in at least one; '
                 'Tier 2 = called in every cohort by at least one method; Tier 3 = called '
                 'in at least two cohorts; 0 = tested in this stratum but called in fewer '
                 'than two cohorts. The tier shown is the stratum\'s own; a gene is listed '
                 'because it reaches Tier 1 or 2 in at least one stratum.\n\n')
        fh.write('**Caveat.** The three cohorts are nested and differ almost entirely in '
                 'controls; agreement across them shows the signal survives adding controls '
                 'and is a sample-selection sensitivity analysis, not replication. Cases and '
                 'controls were sequenced on different platforms, so no gene here is '
                 'established as a CTEPH gene.\n')

    # ---- per-variant and per-sample detail for the tiered genes.
    vs_rows, vs_head = read_tsv(args.variant_stats_all, VARIANT_KEY, 'variant_stats_all.tsv')
    keep = [r for r in vs_rows if (r['stratum'], r['set_name']) in tiered]
    write_table(args.out_variant_detail, vs_head, keep)
    mc_rows, mc_head = read_tsv(args.multi_carriers_all, MULTI_KEY, 'multi_carriers_all.tsv')
    keep_mc = [r for r in mc_rows if (r['stratum'], r['set_name']) in tiered]
    write_table(args.out_multi_detail, mc_head, keep_mc)

    # ---- the windows the per-gene figures draw gene models over.
    windows, seen_w = [], set()
    for sn in order:
        pos = [int(r['pos']) for r in keep
               if r['set_name'] == sn and r['cohort'] == rc and int(r['pos'])]
        if not pos:
            pos = [int(r['pos']) for r in keep if r['set_name'] == sn and int(r['pos'])]
        if not pos or sn in seen_w:
            continue
        seen_w.add(sn)
        lo, hi = min(pos), max(pos)
        pad = max(int(WINDOW_PAD_FRAC * (hi - lo)), WINDOW_PAD_MIN)
        chrom = next(gene_info[(s, sn)][1] for s in strata if (s, sn) in gene_info)
        windows.append({'gene_symbol': summary_sets[sn], 'set_name': sn, 'chrom': chrom,
                        'start': max(lo - pad, 1), 'end': hi + pad})
    write_table(args.out_gene_windows, WINDOW_COLUMNS, windows)

    # ---- report.
    print(f'[robust_genes] evidence: {len(evidence):,} row(s) = '
          f'{sum(len(universe[s]) for s in strata):,} (stratum, gene) x {len(cohorts)} cohort(s) '
          f'x {len(METHODS)} method(s); states '
          + ' '.join(f'{k}={n_state[k]:,}' for k in STATES) + f' -> {args.out_evidence}')
    for t in (1, 2, 3):
        for s in strata:
            g = [r['gene_symbol'] for r in tiers if r['tier'] == t and r['stratum'] == s]
            print(f'    Tier {t}  {s:<18} {len(g):>3}  {", ".join(g) if g else DASH}')
    print(f'[robust_genes] {len(tiers):,} tiered (gene, stratum) row(s) -> {args.out_tiers}; '
          f'{len(order)} gene(s) in the summary tables -> {args.out_summary_prefix}.<cohort>.tsv '
          f'and {args.out_summary_prefix}.md ({rc})')
    print(f'[robust_genes] {len(keep):,} variant row(s) -> {args.out_variant_detail}; '
          f'{len(keep_mc):,} multi-site carrier row(s) -> {args.out_multi_detail}; '
          f'{len(windows)} gene window(s) -> {args.out_gene_windows}')
    flagged = [(r['gene_symbol'], r['stratum'], r['n_multi_site_carriers_case_max'],
                r['n_multi_site_carriers_control_max']) for r in tiers
               if 'multi_site_carriers' in r['note']]
    if flagged:
        print(f'[robust_genes] NOTE: {len(flagged)} tiered (gene, stratum) carry samples with a '
              f'minor allele at >= 2 sites of the same gene (max per cohort, case / control): '
              + '; '.join(f'{g} {s} {a}/{b}' for g, s, a, b in flagged))


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in robust_genes: {e}', file=sys.stderr)
        sys.exit(1)
