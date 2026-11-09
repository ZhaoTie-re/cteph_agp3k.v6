#!/usr/bin/env python3
"""Step 03 — our controls against the two external Japanese panels.

Only CONTROLS are compared (group AGP3K, the population sample); cases are a
disease-selected series whose HLA frequencies may differ for real reasons. All
comparisons are made on IPD-IMGT P groups (ToMMo's naming). Null alleles (suffix N)
are left out on both sides and the rest renormalised: ToMMo uses DRB4*03:01N to mean
"no DRB4 on this haplotype", which our calls record as absence.

A control contributes to a locus only if every copy it carries has a 2-field allele;
the samples left out, and why, are written to excluded.tsv.

Per (panel, locus):
  r                 Pearson correlation of allele frequencies (kept for continuity)
  D                 discordance = half the sum of |f_ours - f_panel| over alleles: the
                    smallest share of chromosomes that would have to change allele to
                    make the two distributions equal
  D_null            D expected from sampling alone: both sides redrawn N_SIM times
                    from their pooled frequencies, at their own chromosome counts
  excess_D          D - mean(D_null). If our controls and the panel share their true
                    frequencies, this is a lower bound on the share of mis-called
                    chromosomes; it cannot see errors that leave frequencies unchanged
  unconfirmed       share of our chromosomes on an allele the panel does not carry,
                    with its expectation under the same sampling model; split into
                    alleles the panel's dictionary could contain (`_assessable`: first
                    named by IPD-IMGT no later than the newest allele the panel reports)
                    and alleles named after it (`_after_vintage`, which cannot be judged)
Per locus (internal, all reported obligate genes, controls):
  homozygosity      observed vs expected under Hardy-Weinberg from our own 2-field
                    allele counts (Levene's exact expectation), with a 95 % Wilson CI
                    on the ratio and a permutation P for excess homozygosity — the
                    signature of allele dropout — and the implied dropout: the share
                    of heterozygotes reported as homozygotes, (obs - exp) / (1 - exp)
  D_shared          D over the alleles both sides carry (renormalised); with
                    D_shared_high and unconfirmed_high, the same on high-confidence calls
Shift direction: for alleles >= 1 point apart from ToMMo at the five 1KG loci, whether
the independent 1KG JPT sample lies on our side or on ToMMo's (shift_direction.tsv).
Panel ceiling: the same D statistics between ToMMo and 1KG JPT on their five shared
loci, i.e. how far two independent panels of the same population disagree.
Confidence check: the unconfirmed share by call confidence (step 02) — whether the
external panel agrees with the read-level grading.

Outputs (results/03.panels/): locus_concordance.tsv, allele_concordance.tsv,
homozygosity.tsv, panel_ceiling.tsv, unconfirmed_alleles.tsv,
confidence_vs_panel.tsv, excluded.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import nomenclature as N  # noqa: E402
from panels import kg_counts, kg_vintage, tommo_counts, tommo_vintage  # noqa: E402

CALLED = ('heterozygous', 'homozygous', 'present')
RNG = np.random.default_rng(C.SEED)


# ------------------------------------------------------------------ inputs
def our_chromosomes(ev, ids):
    """{gene: list of (sample, [P-group alleles])} for the given samples, plus exclusions."""
    out, excl = {}, []
    sub = ev[ev.sample_id.isin(ids)]
    for gene, s in sub.groupby('gene'):
        rows = []
        for r in s.itertuples(index=False):
            if r.state == 'absent':
                continue
            if r.state not in CALLED:
                excl.append({'gene': gene, 'sample_id': r.sample_id, 'reason': r.state})
                continue
            al = [r.allele_1] if r.copies == 1 else [r.allele_1, r.allele_2]
            if any(not a for a in al):
                excl.append({'gene': gene, 'sample_id': r.sample_id, 'reason': 'no_2field'})
                continue
            rows.append((r.sample_id, al, r.confidence))
        out[gene] = rows
    return out, excl


def counts_from(rows):
    """P-group counts (nulls dropped) and the number of null chromosomes dropped."""
    cnt, nulls = {}, 0
    for _, al, _ in rows:
        for a in al:
            if N.is_null(a):
                nulls += 1
                continue
            k = N.pg_of(a)
            cnt[k] = cnt.get(k, 0) + 1
    return cnt, nulls


# ------------------------------------------------------------------ statistics
def disc(a, b):
    na, nb = a.sum(), b.sum()
    return 0.5 * np.abs(a / na - b / nb).sum()


def compare(c_ours, c_panel, n_sim=C.N_SIM):
    alleles = sorted(set(c_ours) | set(c_panel))
    a = np.array([c_ours.get(x, 0) for x in alleles], float)
    b = np.array([c_panel.get(x, 0) for x in alleles], float)
    na, nb = int(a.sum()), int(b.sum())
    d_obs = disc(a, b)
    unconf = a[b == 0].sum() / na
    pooled = (a + b) / (na + nb)
    sa = RNG.multinomial(na, pooled, size=n_sim)
    sb = RNG.multinomial(nb, pooled, size=n_sim)
    d_null = 0.5 * np.abs(sa / na - sb / nb).sum(axis=1)
    u_null = np.array([sa[i][sb[i] == 0].sum() / na for i in range(n_sim)])
    fa, fb = a / na, b / nb
    return {
        'n_chr_ours': na, 'n_chr_panel': nb, 'n_alleles': len(alleles),
        'n_shared': int(((a > 0) & (b > 0)).sum()),
        'r': float(stats.pearsonr(fa, fb)[0]) if len(alleles) > 2 else np.nan,
        'D': d_obs, 'D_null_mean': d_null.mean(), 'D_null_p95': np.percentile(d_null, 95),
        'excess_D': d_obs - d_null.mean(), 'p_D': (1 + (d_null >= d_obs).sum()) / (1 + n_sim),
        'unconfirmed': unconf, 'unconfirmed_null_mean': u_null.mean(),
    }, alleles, a, b


def decompose(c_ours, c_panel):
    """D among the alleles both sides carry, each side renormalised to them."""
    sh = [k for k in c_ours if c_ours[k] and c_panel.get(k, 0)]
    if len(sh) < 2:
        return np.nan
    a = np.array([c_ours[k] for k in sh], float)
    b = np.array([c_panel[k] for k in sh], float)
    return disc(a, b)


def unconfirmed_split(rows, ref, vintage):
    """(chromosomes, on an allele absent from the panel, of which the allele was named by
    the panel's dictionary vintage). Only the last can count against our calls: an allele
    first named after the panel's newest allele cannot appear in it."""
    n = u = ua = 0
    for _, al, _ in rows:
        for a in al:
            if N.is_null(a):
                continue
            n += 1
            if not ref.get(N.pg_of(a), 0):
                u += 1
                r = N.release_of(a)
                if r is not None and r <= vintage:
                    ua += 1
    return n, u, ua


def per_allele(panel, gene, alleles, a, b):
    na, nb = a.sum(), b.sum()
    rows = []
    for x, ka, kb in zip(alleles, a, b):
        p = stats.fisher_exact([[ka, na - ka], [kb, nb - kb]])[1]
        rows.append({'panel': panel, 'gene': gene, 'allele': x, 'count_ours': int(ka),
                     'freq_ours': ka / na, 'count_panel': int(kb), 'freq_panel': kb / nb,
                     'diff': ka / na - kb / nb, 'fisher_p': p})
    df = pd.DataFrame(rows)
    if len(df):
        df['q'] = stats.false_discovery_control(df.fisher_p.values)
    return df


def wilson(k, n, z=1.959963984540054):
    p, d = k / n, 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (c - h) / d, (c + h) / d


def homozygosity(rows, tommo_cnt=None, n_perm=C.N_SIM):
    """Observed vs Hardy-Weinberg homozygosity at the 2-field level (own frequencies)."""
    geno = [(al[0], al[1]) for _, al, _ in rows if len(al) == 2]
    n = len(geno)
    if n < 20:
        return None
    alle = np.array([x for g in geno for x in g])
    obs = sum(1 for x, y in geno if x == y)
    _, c = np.unique(alle, return_counts=True)
    exp = (c * (c - 1)).sum() / (2 * (2 * n - 1))          # Levene: expected homozygotes
    perm = np.empty(n_perm)
    for i in range(n_perm):
        s = RNG.permutation(alle)
        perm[i] = (s[0::2] == s[1::2]).sum()
    lo, hi = wilson(obs, n)
    e = exp / n
    out = {'n_samples': n, 'n_alleles': len(c), 'obs_hom': obs / n, 'exp_hom': e,
           'ratio': obs / exp if exp else np.nan,
           'ratio_lo': lo * n / exp if exp else np.nan, 'ratio_hi': hi * n / exp if exp else np.nan,
           # share of heterozygotes reported as homozygotes, if Hardy-Weinberg holds
           'dropout': (obs / n - e) / (1 - e), 'dropout_lo': (lo - e) / (1 - e),
           'dropout_hi': (hi - e) / (1 - e),
           'p_excess': (1 + (perm >= obs).sum()) / (1 + n_perm)}
    if tommo_cnt:
        tot = sum(tommo_cnt.values())
        out['exp_hom_tommo'] = sum((v / tot) ** 2 for v in tommo_cnt.values())
    return out


# ------------------------------------------------------------------ main
def main():
    ev = pd.read_csv(C.D_EVIDENCE / 'call_evidence.tsv', sep='\t', dtype=str).fillna('')
    ev['copies'] = pd.to_numeric(ev.copies, errors='coerce')
    si = pd.read_csv(C.D_INPUTS / 'sample_inputs.tsv', sep='\t', dtype=str)
    ctrl = set(si.sample_id[si.group == C.CONTROL_GROUP])
    rows, excl = our_chromosomes(ev, ctrl)
    tommo, tommo_nulls = tommo_counts()
    kg, kg_dropped = kg_counts()
    vintage = {'ToMMo': tommo_vintage(), '1KG JPT': kg_vintage()}
    print(f'[panels] dictionary vintage (newest IPD-IMGT release among the panel alleles): {vintage}')

    loc, alle = [], []
    ours_cnt, ours_nulls = {}, {}
    for g in C.REPORTED_GENES:
        ours_cnt[g], ours_nulls[g] = counts_from(rows.get(g, []))
    for panel, ref, loci in (('ToMMo', tommo, C.TOMMO_LOCI), ('1KG JPT', kg, C.KG_LOCI)):
        for g in loci:
            m, al, a, b = compare(ours_cnt[g], ref[g])
            note = C.TOMMO_FREQ_NOT_COMPARABLE.get(g, '') if panel == 'ToMMo' else ''
            hi_rows = [x for x in rows.get(g, []) if x[2] == 'high']
            c_hi, _ = counts_from(hi_rows)
            n_hi = sum(c_hi.values())
            m['D_shared'] = decompose(ours_cnt[g], ref[g])
            m['D_shared_high'] = decompose(c_hi, ref[g])
            m['n_chr_high'] = n_hi
            n_a, u_a, ua_a = unconfirmed_split(rows.get(g, []), ref[g], vintage[panel])
            n_h, u_h, ua_h = unconfirmed_split(hi_rows, ref[g], vintage[panel])
            m['unconfirmed_assessable'] = ua_a / n_a if n_a else np.nan
            m['unconfirmed_after_vintage'] = (u_a - ua_a) / n_a if n_a else np.nan
            m['unconfirmed_high'] = u_h / n_h if n_h else np.nan
            m['unconfirmed_assessable_high'] = ua_h / n_h if n_h else np.nan
            m['panel_vintage'] = vintage[panel]
            loc.append({'panel': panel, 'gene': g, 'n_samples': len(rows.get(g, [])), **m,
                        'null_chr_ours': ours_nulls[g],
                        'null_chr_panel': tommo_nulls.get(g, 0) if panel == 'ToMMo' else '',
                        'panel_values_dropped': kg_dropped.get(g, '') if panel != 'ToMMo' else '',
                        'frequency_comparable': int(not note), 'note': note})
            alle.append(per_allele(panel, g, al, a, b))
    locdf = pd.DataFrame(loc)
    alldf = pd.concat(alle, ignore_index=True)

    # panel ceiling: ToMMo vs 1KG JPT, beside ours vs 1KG JPT
    ceil = []
    for g in C.KG_LOCI:
        m, *_ = compare(tommo[g], kg[g])
        o = locdf[(locdf.panel == '1KG JPT') & (locdf.gene == g)].iloc[0]
        ceil.append({'gene': g, 'D_tommo_vs_kg': m['D'], 'excess_D_tommo_vs_kg': m['excess_D'],
                     'r_tommo_vs_kg': m['r'], 'D_ours_vs_kg': o.D,
                     'excess_D_ours_vs_kg': o.excess_D, 'r_ours_vs_kg': o.r})
    ceildf = pd.DataFrame(ceil)

    # Alleles whose frequency in our controls differs from ToMMo by >= 1 point: is the
    # independent 1KG JPT sample on our side or on ToMMo's?
    t = alldf[alldf.panel == 'ToMMo'][['gene', 'allele', 'freq_ours', 'freq_panel']] \
        .rename(columns={'freq_panel': 'freq_tommo'})
    k = alldf[alldf.panel == '1KG JPT'][['gene', 'allele', 'freq_panel', 'count_panel']] \
        .rename(columns={'freq_panel': 'freq_kg', 'count_panel': 'count_kg'})
    sh = t[t.gene.isin(C.KG_LOCI)].merge(k, on=['gene', 'allele'], how='left').fillna({'freq_kg': 0, 'count_kg': 0})
    sh = sh[(sh.freq_ours - sh.freq_tommo).abs() >= 0.01].copy()
    sh['diff_ours'] = sh.freq_ours - sh.freq_tommo
    sh['diff_kg'] = sh.freq_kg - sh.freq_tommo
    sh['kg_side'] = np.where(np.sign(sh.diff_kg) == np.sign(sh.diff_ours), 'ours', 'ToMMo')

    # homozygosity (internal; obligate reported genes)
    hz = []
    for g in C.OBLIGATE_GENES:
        h = homozygosity(rows.get(g, []), tommo.get(g))
        if h:
            hz.append({'gene': g, **h})
    hzdf = pd.DataFrame(hz)

    # alleles we call that ToMMo does not carry, one row per called 2-field allele, dated
    # by that allele (not by its P group), as in the locus metric
    urows = {}
    for g in C.TOMMO_LOCI:
        for _, al, cf in rows.get(g, []):
            for a in al:
                if N.is_null(a) or tommo[g].get(N.pg_of(a), 0):
                    continue
                r = urows.setdefault((g, a), {'gene': g, 'allele': a, 'p_group': N.pg_label(a),
                                              'count_ours': 0, 'count_high': 0,
                                              'first_release': N.release_of(a)})
                r['count_ours'] += 1
                r['count_high'] += cf == 'high'
    unc = pd.DataFrame(list(urows.values()),
                       columns=['gene', 'allele', 'p_group', 'count_ours', 'count_high', 'first_release'])
    unc['assessable'] = (unc.first_release <= vintage['ToMMo']).astype(int)
    unc = unc.sort_values(['gene', 'count_ours'], ascending=[True, False])

    # unconfirmed share by call confidence (and by the two read flags)
    conf = []
    for g in C.TOMMO_LOCI:
        known = {k for k, v in tommo[g].items() if v > 0}
        for _, al, cf in rows.get(g, []):
            for a in al:
                if N.is_null(a):
                    continue
                u = N.pg_of(a) not in known
                r = N.release_of(a)
                # alleles named after the panel's vintage are left out: they cannot be judged
                if u and not (r is not None and r <= vintage['ToMMo']):
                    continue
                conf.append({'gene': g, 'confidence': cf, 'unconfirmed': u})
    cdf = pd.DataFrame(conf)
    cv = (cdf.groupby('confidence').unconfirmed.agg(['size', 'sum', 'mean'])
          .rename(columns={'size': 'n_chr', 'sum': 'n_unconfirmed', 'mean': 'unconfirmed_share'})
          .reset_index())
    cv.insert(0, 'gene', 'all ToMMo loci')
    cvg = (cdf.groupby(['gene', 'confidence']).unconfirmed.agg(['size', 'sum', 'mean'])
           .rename(columns={'size': 'n_chr', 'sum': 'n_unconfirmed', 'mean': 'unconfirmed_share'})
           .reset_index())

    C.D_PANELS.mkdir(parents=True, exist_ok=True)
    locdf.to_csv(C.D_PANELS / 'locus_concordance.tsv', sep='\t', index=False, float_format='%.6g')
    alldf.to_csv(C.D_PANELS / 'allele_concordance.tsv', sep='\t', index=False, float_format='%.6g')
    hzdf.to_csv(C.D_PANELS / 'homozygosity.tsv', sep='\t', index=False, float_format='%.6g')
    ceildf.to_csv(C.D_PANELS / 'panel_ceiling.tsv', sep='\t', index=False, float_format='%.6g')
    sh.to_csv(C.D_PANELS / 'shift_direction.tsv', sep='\t', index=False, float_format='%.6g')
    unc.to_csv(C.D_PANELS / 'unconfirmed_alleles.tsv', sep='\t', index=False, float_format='%.6g')
    pd.concat([cv, cvg]).to_csv(C.D_PANELS / 'confidence_vs_panel.tsv', sep='\t', index=False,
                                float_format='%.6g')
    pd.DataFrame(excl, columns=['gene', 'sample_id', 'reason']).to_csv(
        C.D_PANELS / 'excluded.tsv', sep='\t', index=False)

    pd.set_option('display.width', 200)
    print(locdf[['panel', 'gene', 'n_samples', 'n_chr_ours', 'r', 'D', 'D_null_mean', 'excess_D',
                 'D_shared', 'D_shared_high', 'unconfirmed', 'unconfirmed_assessable',
                 'unconfirmed_after_vintage', 'unconfirmed_assessable_high']]
          .round(4).to_string(index=False))
    print(ceildf.round(4).to_string(index=False))
    print('[panels] alleles shifted >= 1 point vs ToMMo at the 1KG loci; 1KG JPT on the side of:',
          sh.kg_side.value_counts().to_dict())
    print(hzdf.round(3).to_string(index=False))
    print(cv.round(4).to_string(index=False))


if __name__ == '__main__':
    main()
