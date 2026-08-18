#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Reduce one cohort's three GLM scans to the small tables everything
#           downstream reads.
#
#             scan_qc.tsv        per model: cohort size, lambda_GC, ERRCODE
#                                breakdown, degenerate-fit count, n analysed,
#                                min P.
#             model_peaks.tsv    two-tier peaks for ALL THREE models, merged by
#                                distance. This is what each scan figure draws
#                                itself from, so DOM and REC show their OWN peaks
#                                rather than ADD's.
#             peaks.tsv          ADDITIVE peaks, two tiers, defined by DISTANCE
#                                merging (no LD clumping)
#             lead_variants.tsv  one lead per peak — the fan-out key
#             <peak>.sumstat.tsv per-peak statistics, genome-wide tier only
#
#           Peaks are called ONCE at the suggestive threshold and then tiered,
#           because the two thresholds are nested: every genome-wide peak is a
#           suggestive peak whose strongest variant happens to clear 5e-8.
#           Calling twice would double-count the same physical peak.
#
#           ADD remains the only model that DEFINES the fan-out (nothing downstream
#           reads model_peaks.tsv). Its additive rows therefore duplicate peaks.tsv
#           on purpose: peaks.tsv must stay byte-identical to keep the expensive
#           LD/SuSiE tasks cached, so the general table is written alongside it
#           rather than replacing it.
# Component: assoc_plink2
# Used by : assoc_plink2.nf  process SCAN_PEAKS
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

USECOLS = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'A1_FREQ',
           'OBS_CT', 'OR', 'LOG(OR)_SE', 'L95', 'U95', 'Z_STAT', 'P', 'ERRCODE']
CHI2_MEDIAN = 0.4549364231195729          # qchisq(0.5, df=1)

TIER_GW, TIER_SUG = 'genome_wide', 'suggestive'


def parse_args():
    p = argparse.ArgumentParser(description="Scan QC and additive peak calling for one cohort.")
    p.add_argument('--glm', action='append', required=True, metavar='MODEL=PATH')
    p.add_argument('--cohort', required=True)
    p.add_argument('--fam', help='plink .fam — case/control counts for scan_qc.tsv.')
    p.add_argument('--peak-model', default='additive',
                   help='The ONLY model peaks are called from (default: additive).')
    p.add_argument('--p-genomewide', type=float, default=5e-8)
    p.add_argument('--p-suggestive', type=float, default=1e-5)
    p.add_argument('--peak-flank', type=int, default=250_000,
                   help='Half-width in bp used to merge significant variants into a peak.')
    p.add_argument('--out-dir', required=True)
    return p.parse_args()


def read_glm(path):
    """Read one plink2 GLM file with the C parser (the files are TAB-delimited)."""
    head = pd.read_csv(path, sep='\t', nrows=0)
    cols = [c for c in USECOLS if c in head.columns]
    df = pd.read_csv(path, sep='\t', usecols=cols,
                     dtype={'#CHROM': str, 'ID': str, 'ERRCODE': str})
    df = df.rename(columns={'#CHROM': 'CHROM', 'LOG(OR)_SE': 'SE'})
    for c in ('POS', 'P', 'OR', 'SE', 'L95', 'U95', 'A1_FREQ', 'OBS_CT'):
        if c in df:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df


def usable_mask(df):
    """Rows that carry a result at all. Two conditions, both necessary.

      1. ERRCODE == '.'  — plink2 itself reports a clean fit.
      2. Finite, positive OR and SE, and P strictly inside (0, 1].

    The second is NOT redundant. plink2 leaves ERRCODE '.' on some degenerate
    fits: at a near-fixed variant (A1_FREQ ~0.99) the dominant coding produces
    complete separation and returns OR = inf, SE ~ 4.2e6, |Z| ~ 1e7 and P
    underflowed to exactly 0. Because 0 < 5e-8, such a row satisfies any
    significance test and becomes a phantom genome-wide peak. Three of them did
    exactly that under a non-additive coding before this guard existed.
    """
    err = df['ERRCODE'].fillna('.') if 'ERRCODE' in df else pd.Series('.', index=df.index)
    orv = df['OR'] if 'OR' in df else pd.Series(1.0, index=df.index)
    sev = df['SE'] if 'SE' in df else pd.Series(1.0, index=df.index)
    return (err.eq('.') & df['P'].notna() & (df['P'] > 0) & (df['P'] <= 1)
            & np.isfinite(orv) & (orv > 0) & np.isfinite(sev) & (sev > 0))


def lambda_gc(p):
    """Genomic-control inflation factor from the median chi-square."""
    p = p[(p > 0) & (p <= 1)]
    if not len(p):
        return np.nan
    from scipy import stats
    return float(stats.chi2.isf(np.median(p), df=1) / CHI2_MEDIAN)


def merge_by_distance(sig, flank):
    """Group significant variants into peaks by physical distance.

    Deliberately not LD clumping: a distance rule needs no reference panel, is
    reproducible from the summary statistics alone, and cannot be perturbed by
    the choice of LD sample — which matters here because three different LD
    sources are compared downstream and none of them may define the peak.

    A gap larger than 2*flank starts a new peak; the reported window is the
    variant span extended by `flank` on each side.
    """
    out = []
    for chrom, g in sig.sort_values(['CHROM', 'POS']).groupby('CHROM', sort=False):
        members, start, end = [], None, None
        for _, r in g.iterrows():
            if start is None or r['POS'] - end > 2 * flank:
                if start is not None:
                    out.append((chrom, start, end, members))
                start, end, members = r['POS'], r['POS'], [r]
            else:
                end = r['POS']
                members.append(r)
        if start is not None:
            out.append((chrom, start, end, members))
    return out


def chrom_sort_key(c):
    c = str(c).replace('chr', '')
    return {'X': 23, 'Y': 24, 'MT': 25, 'M': 25}.get(c, int(c) if c.isdigit() else 99)


def cohort_size(fam_path):
    """(n_case, n_ctrl, n_eff) from a plink .fam; phenotype 2 = case, 1 = control.

    N_eff = 4/(1/N_case + 1/N_ctrl) is reported here as a cohort descriptor — the
    balanced-design equivalent size that a later meta-analysis or mixed-model run
    quotes. It used to reach the comparison figure through the power component,
    which no longer exists; scan QC is where a sample size belongs anyway.
    """
    if not fam_path or not Path(fam_path).exists():
        return np.nan, np.nan, np.nan
    ph = pd.read_csv(fam_path, sep=r'\s+', header=None, usecols=[5], names=['PHENO'])
    n1, n0 = int((ph.PHENO == 2).sum()), int((ph.PHENO == 1).sum())
    ne = 4.0 / (1.0 / n1 + 1.0 / n0) if n1 and n0 else np.nan
    return n1, n0, ne


def peaks_by_model(df, alpha_sug, alpha_gw, flank, cohort, model):
    """Two-tier peaks for ONE model, merged into loci by physical distance.

    Identical rule to the additive peak caller — called once at the suggestive
    threshold, then tiered by whether the lead clears 5e-8 — so a peak means the
    same thing in every model and the additive rows here reproduce peaks.tsv
    exactly. Nothing fans out from this table: it exists so every scan's figure can
    be drawn from its OWN result, and so a DOM/REC locus is reportable at all.
    """
    sig = df[usable_mask(df) & (df['P'] < alpha_sug)]
    if not len(sig):
        return []
    groups = merge_by_distance(sig, flank)
    groups.sort(key=lambda t: (chrom_sort_key(t[0]), t[1]))
    rows = []
    for i, (chrom, start, end, members) in enumerate(groups, start=1):
        mem = pd.DataFrame(members)
        lead = mem.loc[mem['P'].idxmin()]
        tier = TIER_GW if lead['P'] < alpha_gw else TIER_SUG
        rows.append({
            'cohort': cohort, 'model': model,
            'peak_id': f"{model[:3]}{i:03d}_{chrom}_{int(lead['POS'])}",
            'tier': tier, 'chrom': chrom,
            'start': int(start - flank), 'end': int(end + flank),
            'pos': int(lead['POS']),
            'n_sig_variants': len(mem),
            'n_genomewide_variants': int((mem['P'] < alpha_gw).sum()),
            'lead_id': lead['ID'], 'lead_p': lead['P'], 'lead_or': lead.get('OR', np.nan),
            'lead_se': lead.get('SE', np.nan),
            'lead_l95': lead.get('L95', np.nan), 'lead_u95': lead.get('U95', np.nan),
            'lead_a1': lead.get('A1', ''), 'lead_ref': lead.get('REF', ''),
            'lead_alt': lead.get('ALT', ''),
            'lead_a1_freq': lead.get('A1_FREQ', np.nan),
            'lead_obs_ct': lead.get('OBS_CT', np.nan),
        })
    return rows


def main():
    args = parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    scans = {}
    for spec in args.glm:
        model, _, path = spec.partition('=')
        if not path:
            raise SystemExit(f'--glm expects MODEL=PATH, got {spec!r}')
        scans[model] = read_glm(path)
    if args.peak_model not in scans:
        raise SystemExit(f'peak model {args.peak_model!r} not among {list(scans)}')

    n_case, n_ctrl, n_eff = cohort_size(args.fam)

    # ── scan QC, all models ─────────────────────────────────────────────────
    qc = []
    for model, df in scans.items():
        err = df['ERRCODE'].fillna('.') if 'ERRCODE' in df else pd.Series('.', index=df.index)
        ok = err.eq('.')
        usable = usable_mask(df)
        clean = df[usable]
        counts = err[~ok].value_counts()
        qc.append({
            'cohort': args.cohort, 'model': model,
            'n_case': n_case, 'n_ctrl': n_ctrl,
            'n_eff': round(n_eff, 1) if np.isfinite(n_eff) else np.nan,
            'n_variants': len(df),
            'n_errcode': int((~ok).sum()),
            'errcode_breakdown': ';'.join(f'{k}={v}' for k, v in counts.items()) or '.',
            'n_degenerate': int((ok & ~usable).sum()),
            'n_analysed': len(clean),
            'obs_ct_median': int(clean['OBS_CT'].median()) if len(clean) and 'OBS_CT' in clean else np.nan,
            'lambda_gc': round(lambda_gc(clean['P']), 5),
            'n_genomewide': int((clean['P'] < args.p_genomewide).sum()),
            'n_suggestive': int((clean['P'] < args.p_suggestive).sum()),
            'min_p': clean['P'].min() if len(clean) else np.nan,
        })
    qc = pd.DataFrame(qc)
    qc['_o'] = qc['model'].map({'additive': 0, 'dominant': 1, 'recessive': 2}).fillna(9)
    qc.sort_values('_o').drop(columns='_o').to_csv(out / 'scan_qc.tsv', sep='\t', index=False)

    # ── peaks for ALL models, both tiers, in the canonical model order ──────
    hit_cols = ['cohort', 'model', 'peak_id', 'tier', 'chrom', 'start', 'end', 'pos',
                'n_sig_variants', 'n_genomewide_variants', 'lead_id',
                'lead_p', 'lead_or', 'lead_se', 'lead_l95', 'lead_u95', 'lead_a1',
                'lead_ref', 'lead_alt', 'lead_a1_freq', 'lead_obs_ct']
    order = {'additive': 0, 'dominant': 1, 'recessive': 2}
    hit_rows = []
    for model in sorted(scans, key=lambda m: order.get(m, 9)):
        hit_rows += peaks_by_model(scans[model], args.p_suggestive, args.p_genomewide,
                                   args.peak_flank, args.cohort, model)
    pd.DataFrame(hit_rows, columns=hit_cols).to_csv(out / 'model_peaks.tsv', sep='\t', index=False)

    # ── peaks: additive only, called once at the suggestive threshold ───────
    add = scans[args.peak_model]
    add = add[usable_mask(add)]
    sig = add[add['P'] < args.p_suggestive]

    peak_rows, lead_rows = [], []
    if len(sig):
        peaks = merge_by_distance(sig, args.peak_flank)
        peaks.sort(key=lambda t: (chrom_sort_key(t[0]), t[1]))
        for i, (chrom, start, end, members) in enumerate(peaks, start=1):
            lo, hi = int(start - args.peak_flank), int(end + args.peak_flank)
            mem = pd.DataFrame(members)
            lead = mem.loc[mem['P'].idxmin()]
            tier = TIER_GW if lead['P'] < args.p_genomewide else TIER_SUG
            prefix = 'gw' if tier == TIER_GW else 'sg'
            peak_id = f"{prefix}{i:03d}_{chrom}_{int(lead['POS'])}"
            ref, alt = lead.get('REF', ''), lead.get('ALT', '')
            row = {
                'cohort': args.cohort, 'peak_id': peak_id, 'tier': tier,
                'chrom': chrom, 'start': lo, 'end': hi,
                'n_sig_variants': len(mem),
                'n_genomewide_variants': int((mem['P'] < args.p_genomewide).sum()),
                'lead_id': lead['ID'], 'lead_pos': int(lead['POS']),
                'lead_p': lead['P'], 'lead_or': lead.get('OR', np.nan),
                'lead_se': lead.get('SE', np.nan),
                'lead_l95': lead.get('L95', np.nan), 'lead_u95': lead.get('U95', np.nan),
                'lead_a1': lead.get('A1', ''), 'lead_ref': ref, 'lead_alt': alt,
                'lead_a1_freq': lead.get('A1_FREQ', np.nan),
                'lead_obs_ct': lead.get('OBS_CT', np.nan),
            }
            peak_rows.append(row)
            lead_rows.append({'cohort': args.cohort, 'peak_id': peak_id, 'tier': tier,
                              'lead_id': lead['ID'], 'chrom': chrom,
                              'pos': int(lead['POS']), 'start': lo, 'end': hi})

    peak_cols = ['cohort', 'peak_id', 'tier', 'chrom', 'start', 'end', 'n_sig_variants',
                 'n_genomewide_variants', 'lead_id', 'lead_pos', 'lead_p', 'lead_or', 'lead_se',
                 'lead_l95', 'lead_u95', 'lead_a1', 'lead_ref', 'lead_alt', 'lead_a1_freq',
                 'lead_obs_ct']
    pd.DataFrame(peak_rows, columns=peak_cols).to_csv(out / 'peaks.tsv', sep='\t', index=False)
    pd.DataFrame(lead_rows, columns=['cohort', 'peak_id', 'tier', 'lead_id', 'chrom',
                                     'pos', 'start', 'end']
                 ).to_csv(out / 'lead_variants.tsv', sep='\t', index=False)

    # ── per-peak summary statistics: genome-wide tier only ──────────────────
    # The suggestive tier gets tables and one display figure, no per-peak
    # follow-up, so writing ~128 window files would be dead weight.
    n_written = 0
    for r in lead_rows:
        if r['tier'] != TIER_GW:
            continue
        win = add[(add['CHROM'] == r['chrom']) & add['POS'].between(r['start'], r['end'])].copy()
        win = win.sort_values('POS')
        win['is_lead'] = win['ID'] == r['lead_id']
        win.to_csv(out / f"{r['peak_id']}.sumstat.tsv", sep='\t', index=False)
        n_written += 1

    n_gw = sum(1 for r in peak_rows if r['tier'] == TIER_GW)
    per_model = {}
    for r in hit_rows:
        k = r['model']
        per_model.setdefault(k, [0, 0])
        per_model[k][0] += 1
        per_model[k][1] += int(r['tier'] == TIER_GW)
    print(f'[call_peaks] {args.cohort}: N={n_case}/{n_ctrl} Neff={n_eff:.1f}; '
          f'{len(peak_rows)} additive peaks ({n_gw} genome-wide, '
          f'{len(peak_rows) - n_gw} suggestive); {n_written} per-peak sumstat file(s); '
          f'per-model peaks ' + (', '.join(f'{m}={n}({g}gw)' for m, (n, g) in per_model.items())
                                 or 'none'))
    print(qc.drop(columns=['errcode_breakdown', '_o'], errors='ignore').to_string(index=False))


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in call_peaks: {e}', file=sys.stderr)
        sys.exit(1)
