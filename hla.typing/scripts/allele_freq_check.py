#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The one external check this component can afford: do our 2-field allele
#           frequencies reproduce the published Japanese ones?
#
#           Everything else in this component's QC measures whether a call was MADE
#           and how deeply it resolved. None of it can see a call that was made
#           confidently and is wrong. A frequency comparison can: if the typing is
#           systematically mis-assigning alleles, the frequency spectrum moves, and
#           it moves in a direction call rate cannot detect.
#
#           The reference is the 1000 Genomes Project's published HLA panel
#           (Abi-Rached et al. 2018), restricted to JPT — the only Japanese
#           population in it, and the reason the file is worth having. It covers
#           five loci at two fields: A, B, C, DQB1, DRB1.
#
#           WHAT THIS IS NOT. It compares POPULATION FREQUENCIES, not genotypes.
#           It cannot say a given sample was typed correctly, and a set of errors
#           that happens to preserve the frequency spectrum is invisible to it. Only
#           typing samples with known types measures accuracy — see
#           docs/OPEN_QUESTIONS.md §2, where that is recorded as deferred.
#
#           TWO CAVEATS ON THE REFERENCE ITSELF:
#             - 105 samples. The standard error on a frequency near 0.4 is ~0.034,
#               so this bounds gross error and nothing finer.
#             - DQB1 is typed on 186 of 210 chromosomes IN THE REFERENCE. Its
#               denominator is reported per locus rather than assumed to be 2n.
#
#           The controls are what is compared. They are an unselected Japanese
#           population sample; the cases are a disease series and are not expected
#           to match a population reference at every locus, least of all in the MHC.
#           Their frequencies are reported beside, never as the comparison.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import re
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

USABLE = ('called', 'hemizygous')

# The reference is NOT uniformly `NN:NN`, whatever its documentation suggests. Of
# its 453 distinct values, 39 carry a trailing '*', 63 are '/'-separated ambiguity
# strings, and a handful are neither an allele nor anything else ('None', a
# comma-decimal number that leaked into an allele column, a space-separated pair).
# Reading them literally invents allele names that can never match a call, which
# both inflates "present in the reference only" and puts junk rows in the table.
REF_CLEAN = re.compile(r'^(\d+:\d+[A-Z]?)\*?$')   # 02:07* -> 02:07, 04:09N kept


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--calls', required=True, help='allele_calls.tsv')
    p.add_argument('--manifest', required=True, help='sample_manifest.tsv, for `group`')
    p.add_argument('--truth', required=True,
                   help='20181129_HLA_types_full_1000_Genomes_Project_panel.txt')
    p.add_argument('--population', default='JPT',
                   help='the reference population. JPT is the only Japanese one.')
    p.add_argument('--control-group', default='AGP3K',
                   help='the `group` value that is the population sample')
    p.add_argument('--genes', default='A,B,C,DQB1,DRB1',
                   help='the five loci the reference carries; nothing else is comparable')
    p.add_argument('--field-depth', type=int, default=2,
                   help="the reference's resolution, and therefore the comparison's")
    p.add_argument('--out', default='allele_frequency_check.tsv')
    p.add_argument('--out-summary', default='allele_frequency_summary.tsv')
    return p.parse_args()


def trim(allele, depth):
    f = str(allele).split(':')
    return ':'.join(f[:depth]) if len(f) > depth else str(allele)


def wilson(k, n, z=1.959963984540054):
    """95 % Wilson score interval. Normal-approximation intervals put the bound
    below 0 for the rare alleles that make up most of an HLA frequency table."""
    if not n:
        return float('nan'), float('nan')
    p, d = k / n, 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5)
    return max(0.0, (c - h) / d), min(1.0, (c + h) / d)


def observed_counts(calls, genes, depth):
    """{gene: (Counter-as-dict, n_chromosomes)} from one allele_calls.tsv slice."""
    out = {}
    for g in genes:
        c1, c2, cs = f'{g}_1', f'{g}_2', f'{g}_state'
        if c1 not in calls.columns:
            out[g] = ({}, 0)
            continue
        counts, n = {}, 0
        for _, r in calls.iterrows():
            if r.get(cs, '') not in USABLE:
                continue
            for a in (r[c1], r[c2]):
                if not a:
                    continue
                # Hemizygous is written as the same allele twice, which is what a
                # homozygote is at these five loci. The DRB3/4/5 defect of
                # OPEN_QUESTIONS §1 does not reach here: none of them is compared.
                counts[trim(a, depth)] = counts.get(trim(a, depth), 0) + 1
                n += 1
        out[g] = (counts, n)
    return out


def reference_counts(truth_path, population, genes, depth):
    """{gene: (counts, n_chromosomes)} for one 1KG population."""
    d = pd.read_csv(truth_path, sep='\t', dtype=str)
    for col in ('Population', 'Sample ID'):
        if col not in d.columns:
            raise SystemExit(f'ABORT: {truth_path} has no {col!r} column. '
                             f'Found: {list(d.columns)[:6]} …')
    sub = d[d['Population'] == population]
    if not len(sub):
        raise SystemExit(f'ABORT: no sample with Population == {population!r} in '
                         f'{truth_path}. Present: '
                         f'{sorted(d["Population"].dropna().unique())}')
    out, dropped = {}, {}
    for g in genes:
        c1, c2 = f'HLA-{g} 1', f'HLA-{g} 2'
        if c1 not in d.columns:
            # The panel carries five loci. A gene it does not carry is recorded as
            # having no reference, never silently compared against nothing.
            out[g] = ({}, 0)
            continue
        counts, n, bad = {}, 0, []
        for a in pd.concat([sub[c1], sub[c2]]):
            if pd.isna(a) or str(a).strip() in ('', '-', 'NA'):
                continue                       # untyped: out of the denominator
            m = REF_CLEAN.match(str(a).strip())
            if not m:
                # An ambiguity string ('35:01/40N/42/57/94') or junk. It cannot be
                # attributed to one allele, so it leaves the DENOMINATOR too — a
                # chromosome counted but never assignable would deflate every
                # frequency at this locus by a little.
                bad.append(str(a).strip())
                continue
            # The panel writes `23:01` with no gene prefix; a trailing '*' marks a
            # call the source flagged, not a different allele, so it is merged in.
            key = f'{g}*{trim(m.group(1), depth)}'
            counts[key] = counts.get(key, 0) + 1
            n += 1
        out[g] = (counts, n)
        dropped[g] = bad
    return out, dropped


def main():
    args = parse_args()
    genes = [g.strip() for g in args.genes.split(',') if g.strip()]

    calls = pd.read_csv(args.calls, sep='\t', dtype=str).fillna('')
    man = pd.read_csv(args.manifest, sep='\t', dtype=str)
    if 'group' not in man.columns:
        raise SystemExit(f'ABORT: {args.manifest} has no `group` column')
    grp = dict(zip(man.sample_id, man.group))
    calls['__group'] = calls.sample_id.map(grp)

    ctrl = calls[calls.__group == args.control_group]
    case = calls[calls.__group != args.control_group]
    # No controls is not an error: a pilot spread over platforms draws its first
    # samples from the case platforms, and there are none. The reference columns are
    # still written, so the table's shape does not depend on the sample set, and
    # n_chr_ctrl == 0 says plainly that nothing was compared. That a FULL run has a
    # usable `group` column is BUILD_MANIFEST's job, and it aborts if the column is
    # missing — that check does not belong here.
    if not len(ctrl):
        print(f'[allele_freq_check] WARNING: 0 sample(s) in group '
              f'{args.control_group!r} (present: '
              f'{sorted(set(calls.__group.dropna())) or "none"}). '
              f'Nothing to compare; writing the reference columns only.',
              file=sys.stderr)

    o_ctrl = observed_counts(ctrl, genes, args.field_depth)
    o_case = observed_counts(case, genes, args.field_depth)
    o_ref, ref_dropped = reference_counts(args.truth, args.population, genes,
                                          args.field_depth)
    n_dropped = sum(len(v) for v in ref_dropped.values())
    if n_dropped:
        print(f'[allele_freq_check] {n_dropped} reference value(s) were not a 2-field '
              f'allele and left both numerator and denominator:', file=sys.stderr)
        for g, vals in ref_dropped.items():
            if vals:
                print(f'    {g}: {", ".join(sorted(set(vals)))}', file=sys.stderr)

    rows, summary = [], []
    for g in genes:
        cc, nc = o_ctrl[g]
        ca, na = o_case[g]
        cr, nr = o_ref[g]
        alleles = sorted(set(cc) | set(cr))
        for al in alleles:
            k_c, k_r = cc.get(al, 0), cr.get(al, 0)
            f_c = k_c / nc if nc else float('nan')
            f_r = k_r / nr if nr else float('nan')
            lo, hi = wilson(k_c, nc)
            # Fisher against the reference COUNTS, not against its frequency: the
            # reference is itself an estimate from 105 samples and treating it as
            # known would make every difference look significant.
            p = (stats.fisher_exact([[k_c, nc - k_c], [k_r, nr - k_r]])[1]
                 if nc and nr else float('nan'))
            rows.append({
                'gene': g, 'allele': al,
                'count_ctrl': k_c, 'n_chr_ctrl': nc,
                'freq_ctrl': round(f_c, 5) if nc else '',
                'ci_lo_ctrl': round(lo, 5) if nc else '',
                'ci_hi_ctrl': round(hi, 5) if nc else '',
                'count_case': ca.get(al, 0), 'n_chr_case': na,
                'freq_case': round(ca.get(al, 0) / na, 5) if na else '',
                'count_ref': k_r, 'n_chr_ref': nr,
                'freq_ref': round(f_r, 5) if nr else '',
                'diff_ctrl_ref': round(f_c - f_r, 5) if nc and nr else '',
                'fisher_p_ctrl_ref': f'{p:.3g}' if p == p else '',
                'in_ref_only': int(k_c == 0 and k_r > 0),
                'in_obs_only': int(k_r == 0 and k_c > 0),
            })
        if nc and nr and len(alleles) > 2:
            fc = [cc.get(a, 0) / nc for a in alleles]
            fr = [cr.get(a, 0) / nr for a in alleles]
            worst = max(range(len(alleles)), key=lambda i: abs(fc[i] - fr[i]))
            summary.append({
                'gene': g, 'n_alleles': len(alleles),
                'n_chr_ctrl': nc, 'n_chr_ref': nr,
                'pearson_r': round(float(stats.pearsonr(fc, fr)[0]), 4),
                'spearman_rho': round(float(stats.spearmanr(fc, fr)[0]), 4),
                'max_abs_diff': round(abs(fc[worst] - fr[worst]), 4),
                'allele_at_max': alleles[worst],
                'n_ref_only': sum(1 for a in alleles if cc.get(a, 0) == 0),
                'n_obs_only': sum(1 for a in alleles if cr.get(a, 0) == 0),
            })
        else:
            summary.append({'gene': g, 'n_alleles': len(alleles),
                            'n_chr_ctrl': nc, 'n_chr_ref': nr,
                            'pearson_r': '', 'spearman_rho': '', 'max_abs_diff': '',
                            'allele_at_max': '', 'n_ref_only': '', 'n_obs_only': ''})

    cols = ['gene', 'allele', 'count_ctrl', 'n_chr_ctrl', 'freq_ctrl', 'ci_lo_ctrl',
            'ci_hi_ctrl', 'count_case', 'n_chr_case', 'freq_case', 'count_ref',
            'n_chr_ref', 'freq_ref', 'diff_ctrl_ref', 'fisher_p_ctrl_ref',
            'in_ref_only', 'in_obs_only']
    pd.DataFrame(rows, columns=cols).to_csv(args.out, sep='\t', index=False)
    scols = ['gene', 'n_alleles', 'n_chr_ctrl', 'n_chr_ref', 'pearson_r',
             'spearman_rho', 'max_abs_diff', 'allele_at_max', 'n_ref_only',
             'n_obs_only']
    pd.DataFrame(summary, columns=scols).to_csv(args.out_summary, sep='\t', index=False)

    print(f'[allele_freq_check] {len(ctrl)} control(s) vs {args.population} '
          f'({o_ref[genes[0]][1] // 2} samples) over {len(genes)} locus/loci')
    for s in summary:
        if s['pearson_r'] == '':
            why = ('no control sample' if not s['n_chr_ctrl'] else
                   'no reference for this locus' if not s['n_chr_ref'] else
                   'too few alleles to correlate')
            print(f"      {s['gene']:6s} n_chr {s['n_chr_ctrl']:>6} vs {s['n_chr_ref']:>4}"
                  f"  — {why}")
        else:
            print(f"      {s['gene']:6s} n_chr {s['n_chr_ctrl']:>6} vs {s['n_chr_ref']:>4}"
                  f"  r={s['pearson_r']:.3f}  rho={s['spearman_rho']:.3f}"
                  f"  max|diff|={s['max_abs_diff']:.3f} at {s['allele_at_max']}")
    print(f'    {len(rows)} (gene, allele) row(s) -> {args.out}')
    print('    This compares FREQUENCIES. It does not measure per-sample accuracy.')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in allele_freq_check: {e}', file=sys.stderr)
        sys.exit(1)
