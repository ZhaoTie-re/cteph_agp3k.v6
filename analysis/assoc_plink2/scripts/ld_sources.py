#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Assemble linkage disequilibrium with the lead variant of one locus
#           from three independent sources, on one common scale:
#             cohort      in-sample, the samples the statistics came from
#             tommo       population co-occurrence r2 tables (tabix)
#             1000g_eas   an out-of-sample reference genotype panel
#           The three keys are stable identifiers for the three ROLES; which
#           resource fills each is configuration (--pop-r2-template,
#           --ref-panel-bfile) and its display name is --panel-label.
#           Also writes the SIGNED square r matrix over the locus variants,
#           which is what SuSiE needs (r2 loses the sign and would make every
#           pair of variants look positively correlated).
#
#           Truncation. Co-occurrence tables typically record pairs with
#           r2 >= 0.2 only. A pair absent from the table is therefore BOUNDED
#           pairs, minimum r2 = 0.200054, none below. Absence of a pair is
#           therefore a BOUND, not missing data. Two states must be separated:
#             - variant is in the panel (its self-pair is present) but has no
#               pair with the lead  ->  r2 < 0.2, the lowest bin, like any other
#               low-LD variant;
#             - variant is not in the panel at all  ->  genuinely unknown, grey.
#           Panel membership is read off the self-pairs the table already
#           carries, so it costs no extra query.
# Component: assoc_plink2
# Used by : assoc_plink2.nf  process LD_SOURCES
# ---------------------------------------------------------------------------
import argparse
import gzip
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MEASURED, BELOW, UNKNOWN = 'measured', 'below_threshold', 'not_in_panel'
POP_R2_FLOOR = 0.2      # co-occurrence tables publish pairs at or above this


def parse_args():
    p = argparse.ArgumentParser(description='LD with the locus lead from three sources.')
    p.add_argument('--plink', required=True, help='PLINK 1.9 — provides --r/--r2 square.')
    p.add_argument('--plink2', required=True, help='PLINK 2 — used to cut windows out of large panels.')
    p.add_argument('--tabix', required=True)
    p.add_argument('--bfile', required=True, help='Cohort fixed_model bfile (in-sample LD).')
    p.add_argument('--ref-panel-bfile', dest='eas_bfile',
                   help='Out-of-sample reference panel bfile prefix; omit to skip that source.')
    p.add_argument('--panel-label', action='append', default=[], metavar='KEY=LABEL',
                   help='Display name per LD source key, recorded in the coverage table.')
    p.add_argument('--pop-r2-template', dest='tommo_template',
                   help="Per-chromosome co-occurrence r2 tables, e.g. '/path/...-chr@@CHROM@@-r2.tsv.gz'")
    p.add_argument('--sumstat', required=True, help='<locus>.sumstat.tsv from SCAN_SUMMARY.')
    p.add_argument('--locus-id', required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--lead-id', required=True)
    p.add_argument('--chrom', required=True)
    p.add_argument('--start', type=int, required=True)
    p.add_argument('--end', type=int, required=True)
    p.add_argument('--memory-mb', type=int, default=8000)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--out-dir', required=True)
    return p.parse_args()


def sh(cmd):
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(' '.join(map(str, cmd)), file=sys.stderr)
        print(res.stdout[-3000:] + res.stderr[-3000:], file=sys.stderr)
    return res.returncode == 0


# ── source 1: in-sample cohort LD ──────────────────────────────────────────
def cohort_ld(args, ids, work):
    """Signed square r matrix over `ids`. --keep-allele-order fixes the sign
    convention to REF/ALT so the matrix and the GLM betas refer to the same
    allele; without it PLINK recodes by minor allele and the signs drift."""
    ex = work / 'locus.vars'
    ex.write_text('\n'.join(ids) + '\n')
    pfx = work / 'cohort_r'
    ok = sh([args.plink, '--bfile', args.bfile, '--extract', str(ex),
             '--keep-allele-order', '--r', 'square', '--memory', str(args.memory_mb),
             '--threads', str(args.threads), '--out', str(pfx)])
    if not ok or not (Path(str(pfx) + '.ld')).exists():
        return None, None
    # PLINK emits the matrix in .bim order over the extracted subset; recover
    # that order from the run's own .nosex/.log-free artefacts by re-reading the
    # bim and filtering, which is the only order PLINK guarantees.
    bim = pd.read_csv(args.bfile + '.bim', sep='\t', header=None,
                      names=['CHROM', 'ID', 'CM', 'POS', 'A1', 'A2'], dtype={'ID': str})
    order = [i for i in bim['ID'] if i in set(ids)]
    mat = np.loadtxt(str(pfx) + '.ld')
    mat = np.atleast_2d(mat)
    if mat.shape[0] != len(order):
        print(f'[ld_sources] cohort matrix {mat.shape} vs {len(order)} ids — abort', file=sys.stderr)
        return None, None
    return pd.DataFrame(mat, index=order, columns=order), order


# ── source 2: population co-occurrence r2 tables ───────────────────────────
def tommo_ld(args, ids, lead_key):
    """Return {variant_key: r2_with_lead} plus the set of keys present in the panel.

    Panel membership comes from the self-pairs (variation1 == variation2) the
    table carries for every panel variant in the queried interval.
    """
    if not args.tommo_template:
        return {}, set()
    path = args.tommo_template.replace('@@CHROM@@', str(args.chrom).replace('chr', ''))
    if not Path(path).exists():
        print(f'[ld_sources] co-occurrence table absent: {path}', file=sys.stderr)
        return {}, set()
    region = f"chr{str(args.chrom).replace('chr', '')}:{args.start}-{args.end}"
    proc = subprocess.Popen([args.tabix, path, region], stdout=subprocess.PIPE, text=True)

    r2, panel = {}, set()
    for line in proc.stdout:
        f = line.rstrip('\n').split('\t')
        if len(f) < 13:
            continue
        k1 = f'chr{f[1].replace("chr", "")}:{f[2]}:{f[3]}:{f[4]}'
        k2 = f'chr{f[7].replace("chr", "")}:{f[8]}:{f[9]}:{f[10]}'
        panel.add(k1)
        panel.add(k2)
        if k1 == lead_key and k2 != lead_key:
            r2[k2] = float(f[12])
        elif k2 == lead_key and k1 != lead_key:
            r2[k1] = float(f[12])
    proc.wait()
    # The lead's own self-pair is r2 = 1 by definition.
    if lead_key in panel:
        r2[lead_key] = 1.0
    return r2, panel


# ── source 3: out-of-sample reference genotype panel ───────────────────────
def eas_ld(args, work):
    """r2 against the lead in the out-of-sample reference panel.

    Two steps on purpose. The panel holds 67 M variants, and running --r2 against
    it directly makes PLINK size its workspace from the whole file — a 34 GB
    allocation that fails. Cutting the window out with --make-bed first reduces
    that to a few thousand variants. --ld-snp then yields one row per partner, so
    no square matrix is materialised (the reference project wrote a 40-87 MB
    matrix per locus and never removed it).
    """
    if not args.eas_bfile:
        return {}, set()
    chrom = str(args.chrom).replace('chr', '')
    sub = work / 'eas_sub'
    if not sh([args.plink2, '--bfile', args.eas_bfile, '--chr', chrom,
               '--from-bp', str(args.start), '--to-bp', str(args.end),
               '--make-bed', '--memory', str(args.memory_mb),
               '--threads', str(args.threads), '--out', str(sub)]) \
            or not Path(str(sub) + '.bed').exists():
        return {}, set()
    win_kb = int((args.end - args.start) / 1000) + 10
    n_var = sum(1 for _ in open(str(sub) + '.bim'))
    pfx = work / 'eas_r2'
    ok = sh([args.plink, '--bfile', str(sub),
             '--ld-snp', args.lead_id, '--r2',
             '--ld-window-kb', str(win_kb), '--ld-window', str(n_var + 1),
             '--ld-window-r2', '0',
             '--memory', str(args.memory_mb), '--threads', str(args.threads),
             '--out', str(pfx)])
    # Panel membership from the .bim regardless of whether the lead was found:
    # "the lead is absent" and "the partner is absent" are different statements.
    bim = pd.read_csv(str(sub) + '.bim', sep='\t', header=None,
                      names=['CHROM', 'ID', 'CM', 'POS', 'A1', 'A2'],
                      dtype={'ID': str, 'CHROM': str})
    panel = set(bim['ID'])
    if not ok or not Path(str(pfx) + '.ld').exists():
        return {}, panel
    ld = pd.read_csv(str(pfx) + '.ld', sep=r'\s+')
    out = {}
    for _, r in ld.iterrows():
        other = r['SNP_B'] if r['SNP_A'] == args.lead_id else r['SNP_A']
        out[other] = float(r['R2'])
    out[args.lead_id] = 1.0
    return out, panel


def classify(ids, r2map, panel, floor_is_bound):
    """Map every locus variant to (r2, state) under one source's semantics."""
    rows = []
    for vid in ids:
        if vid in r2map:
            rows.append((vid, r2map[vid], MEASURED))
        elif vid in panel:
            # In the panel but no pair reported. Under the table's reporting floor
            # that is the bounded statement r2 < 0.2 (flagged BELOW so the
            # figure can say so); for a source with no floor it is a true zero.
            rows.append((vid, 0.0, BELOW if floor_is_bound else MEASURED))
        else:
            rows.append((vid, np.nan, UNKNOWN))
    return pd.DataFrame(rows, columns=['ID', 'r2', 'state'])


def main():
    args = parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    work = out / '_work'
    work.mkdir(exist_ok=True)

    ss = pd.read_csv(args.sumstat, sep='\t', dtype={'CHROM': str, 'ID': str})
    ids = ss['ID'].tolist()

    # ── cohort: square signed r (for SuSiE) and r2 vs lead (for plots) ──────
    mat, order = cohort_ld(args, ids, work)
    coverage = []
    if mat is not None:
        np.savetxt(out / f'{args.locus_id}.ld_matrix.tsv',
                   mat.to_numpy(), fmt='%.6f', delimiter='\t')
        pd.Series(order, name='ID').to_csv(out / f'{args.locus_id}.ld_matrix.vars',
                                           index=False, header=True)
        if args.lead_id in mat.index:
            r = mat.loc[args.lead_id]
            coh = pd.DataFrame({'ID': r.index, 'r2': (r.values ** 2), 'state': MEASURED})
        else:
            coh = pd.DataFrame({'ID': ids, 'r2': np.nan, 'state': UNKNOWN})
        coh = coh.set_index('ID').reindex(ids).reset_index()
        coh['state'] = coh['state'].fillna(UNKNOWN)
    else:
        coh = pd.DataFrame({'ID': ids, 'r2': np.nan, 'state': UNKNOWN})
    coh.to_csv(out / f'{args.locus_id}.ld_cohort.tsv', sep='\t', index=False, float_format='%.6f')

    # ── population co-occurrence panel ─────────────────────────────────────
    lead_row = ss[ss['ID'] == args.lead_id]
    lead_key = args.lead_id
    tm_r2, tm_panel = tommo_ld(args, ids, lead_key)
    tm = classify(ids, tm_r2, tm_panel, floor_is_bound=True)
    tm.to_csv(out / f'{args.locus_id}.ld_tommo.tsv', sep='\t', index=False, float_format='%.6f')

    # ── reference panel ────────────────────────────────────────────────────
    eas_r2, eas_panel = eas_ld(args, work)
    eas = classify(ids, eas_r2, eas_panel, floor_is_bound=False)
    eas.to_csv(out / f'{args.locus_id}.ld_1000g_eas.tsv', sep='\t', index=False, float_format='%.6f')

    labels = dict(s.split('=', 1) for s in args.panel_label if '=' in s)
    for name, tbl, note in [
            ('cohort', coh, 'in-sample, all analysed samples'),
            ('tommo', tm, f'{labels.get("tommo", "population panel")} co-occurrence, '
                          f'reporting floor r2 >= {POP_R2_FLOOR}'),
            ('1000g_eas', eas, f'{labels.get("1000g_eas", "reference panel")}, out-of-sample')]:
        c = tbl['state'].value_counts()
        coverage.append({'cohort': args.cohort, 'locus_id': args.locus_id, 'source': name,
                         'n_variants': len(tbl), 'n_measured': int(c.get(MEASURED, 0)),
                         'n_below_threshold': int(c.get(BELOW, 0)),
                         'n_not_in_panel': int(c.get(UNKNOWN, 0)),
                         'pct_informative': round(100 * (len(tbl) - c.get(UNKNOWN, 0)) / max(len(tbl), 1), 1),
                         'note': note})
    pd.DataFrame(coverage).to_csv(out / f'{args.locus_id}.ld_coverage.tsv', sep='\t', index=False)

    for pat in ('*.ld', '*.nosex', 'eas_sub.*'):
        for f in work.glob(pat):
            f.unlink()
    print(f'[ld_sources] {args.cohort} {args.locus_id}: {len(ids)} variants')
    print(pd.DataFrame(coverage).drop(columns=['note']).to_string(index=False))


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in ld_sources: {e}', file=sys.stderr)
        sys.exit(1)
