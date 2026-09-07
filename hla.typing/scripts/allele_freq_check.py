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
#           The reference is ToMMo's jMorp 61KJPN-HLA panel: 61,424 Japanese
#           individuals over 13 loci. The 1000 Genomes JPT panel it replaced is
#           still supported (--truth-format 1kg_wide) and kept as an independent
#           second reference; the two agree at r = 0.94-0.98 on the five loci they
#           share, which is worth more than either alone.
#
#           WHY IT WAS REPLACED. JPT is 105 individuals over A, B, C, DQB1 and
#           DRB1. At n = 105 the standard error near a frequency of 0.4 is ~0.034,
#           so it bounds gross error and nothing finer — and it carries no DPB1,
#           no DQA1 and no DRB3/4/5, which are the loci this pipeline types least
#           reliably. The new panel is ~580x larger and covers 13 loci.
#
#           P GROUPS. jMorp names alleles by IPD-IMGT P group (A*01:01P) and
#           HLA-HD does not. Same alleles, two spellings; joining without the
#           translation is wrong rather than coarse, because DRB4*01:03 — 76 % of
#           our DRB4 chromosomes — is a member of DRB4*01:01P and untranslated
#           reads as absent from a 61,424-person panel. The translation is applied
#           to BOTH sides here and NOWHERE ELSE: allele_dosage.tsv and
#           residue_dosage.tsv keep plain 2-field names, because P groups merge
#           alleles differing outside the antigen recognition domain and the
#           residue analysis reads the full protein.
#
#           WHAT THIS IS NOT. It compares POPULATION FREQUENCIES, not genotypes.
#           It cannot say a given sample was typed correctly, and a set of errors
#           that happens to preserve the frequency spectrum is invisible to it. Only
#           typing samples with known types measures accuracy — see
#           docs/OPEN_QUESTIONS.md §2, where that is recorded as deferred.
#
#           CAVEATS ON THE REFERENCE ITSELF:
#             - Denominators are read per locus, never assumed to be 2n. jMorp's
#               are uniform; 1KG's are not (its DQB1 is typed on 186 of 210).
#             - DRB4 encodes gene absence as the null allele DRB4*03:01N, 0.5430
#               of all chromosomes, which is what this pipeline calls `not_typed`.
#               Null alleles are dropped and the rest renormalised, or the two
#               sides are not the same quantity. That is the whole DRB4 gap:
#               dropping them takes it from r = 0.625 to r = 0.991.
#             - DRB3 does NOT reconcile and is not made to. It has no null alleles
#               yet spans all 122,848 chromosomes, and sits at r ~ 0.64. The
#               translation, null removal and the OPEN_QUESTIONS §1 defect were all
#               tested and none explains it. Reported, not guessed at.
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
                   help='the reference panel; see --truth-format')
    p.add_argument('--truth-format', default='jmorp_long',
                   choices=('jmorp_long', '1kg_wide'),
                   help="'jmorp_long' = jMorp's gene,allele,count,frequency table over "
                        "61,424 Japanese. '1kg_wide' = the 1000 Genomes per-sample panel, "
                        "kept as an independent second reference over its five loci.")
    p.add_argument('--p-groups', default=None,
                   help="IPD-IMGT/HLA hla_nom_p.txt. REQUIRED for jmorp_long: jMorp names "
                        "alleles by P group and we do not, so without it our commonest "
                        "DRB4 call reads as absent from a 61,424-person panel. Applied to "
                        "BOTH sides so the comparison is like for like, and NEVER to the "
                        "dosage tables the association reads.")
    p.add_argument('--population', default='JPT',
                   help='1kg_wide only: the reference population. JPT is the only Japanese one.')
    p.add_argument('--control-group', default='AGP3K',
                   help='the `group` value that is the population sample')
    p.add_argument('--genes', default='A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1,DRB3,DRB4,E,F,G',
                   help='the loci the reference carries; nothing else is comparable')
    p.add_argument('--field-depth', type=int, default=2,
                   help="the reference's resolution, and therefore the comparison's")
    p.add_argument('--out', default='allele_frequency_check.tsv')
    p.add_argument('--out-pgroup-map', default='allele_pgroup_map.tsv',
                   help='crosswalk for the association that follows: every allele this '
                        'cohort actually carries, its P group, and the reference frequency '
                        'behind it. Emitted HERE rather than from residue_matrix.py so '
                        'that allele_dosage.tsv and residue_dosage.tsv are not '
                        'regenerated — the association reads those, and they must not '
                        'move because a QC reference changed.')
    p.add_argument('--out-summary', default='allele_frequency_summary.tsv')
    p.add_argument('--out-sample', default='allele_frequency_sample.tsv',
                   help='per-sample counts of chromosomes on a reference-confirmed allele. '
                        'This is the only per-SAMPLE accuracy proxy the component has: call '
                        'rate and field depth are near-saturated and say nothing about whether '
                        'a call is right.')
    return p.parse_args()


def trim(allele, depth):
    f = str(allele).split(':')
    return ':'.join(f[:depth]) if len(f) > depth else str(allele)


def load_pgroups(path):
    """{'GENE*NN:NN': 'GENE*NN:NNP'} from IPD-IMGT/HLA's hla_nom_p.txt.

    A P group collects the alleles whose protein sequence is identical over the
    ANTIGEN RECOGNITION DOMAIN — exons 2+3 for class I, exon 2 for class II — and is
    named after its lowest member. jMorp reports frequencies against these names;
    HLA-HD reports plain allele names. They are the same alleles spelled two ways,
    and joining without the translation is wrong rather than imprecise: DRB4*01:03 is
    76 % of our DRB4 chromosomes and is a member of DRB4*01:01P, so untranslated it
    reads as absent from a 61,424-person panel.

    File format, one group per line:  GENE*;a1/a2/.../aN;PGROUP
    PGROUP is empty when the allele is alone in its group.

    Keyed at 2 fields because that is the depth the comparison runs at, which is only
    legitimate if a 2-field name never spans two P groups. On 3.64.0 none of the
    25,290 of them does, and that is asserted below rather than assumed.
    """
    two = {}
    conflict = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('#') or ';' not in line:
            continue
        parts = line.split(';')
        if len(parts) < 3:
            continue
        gene, members, pg = parts[0].rstrip('*'), parts[1], parts[2].strip()
        for a in members.split('/'):
            key = f"{gene}*{':'.join(a.split(':')[:2])}"
            val = f'{gene}*{pg}' if pg else key
            if two.setdefault(key, val) != val:
                conflict.append((key, two[key], val))
    if conflict:
        # The whole 2-field comparison rests on this being single-valued. If a
        # release ever breaks it the comparison is invalid, so this aborts rather
        # than warning and carrying on with a silently arbitrary choice.
        s = '; '.join(f'{k} -> {a} and {b}' for k, a, b in conflict[:5])
        raise SystemExit(f'ABORT: {len(conflict)} two-field allele name(s) span more than '
                         f'one P group in {path}, so a 2-field comparison is not '
                         f'well defined. First: {s}')
    print(f'[allele_freq_check] P groups: {len(two):,} two-field names, all single-valued')
    return two


def to_pgroup(key, pmap):
    """A 'GENE*NN:NN' key translated to its P group, or left alone if it has none."""
    return pmap.get(key, key) if pmap else key


def wilson(k, n, z=1.959963984540054):
    """95 % Wilson score interval. Normal-approximation intervals put the bound
    below 0 for the rare alleles that make up most of an HLA frequency table."""
    if not n:
        return float('nan'), float('nan')
    p, d = k / n, 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5)
    return max(0.0, (c - h) / d), min(1.0, (c + h) / d)


def observed_counts(calls, genes, depth, pmap=None):
    """{gene: (Counter-as-dict, n_chromosomes)} from one allele_calls.tsv slice.

    `pmap` translates each key to its P group so the counts are on the reference's
    own nomenclature. This is LOCAL to the frequency check: allele_dosage.tsv and
    residue_dosage.tsv keep the untranslated 2-field names, because P groups merge
    alleles that differ outside the antigen recognition domain and the residue
    analysis reads the full protein.
    """
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
                # homozygote is at the classical loci. DRB3/DRB4 ARE compared now, so
                # the OPEN_QUESTIONS §1 defect does reach here — measured, it moves the
                # per-locus r by less than 0.02 everywhere (DRB4 0.9905 vs 0.9785,
                # DRB3 0.6417 vs 0.6343). It distorts DOSAGES, not frequencies.
                k = to_pgroup(trim(a, depth), pmap)
                counts[k] = counts.get(k, 0) + 1
                n += 1
        out[g] = (counts, n)
    return out


def reference_counts_jmorp(truth_path, genes, pmap):
    """{gene: (counts, n_chromosomes)} from jMorp's long allele-frequency table.

    Schema `gene,allele,count,frequency` over 61,424 Japanese individuals; every gene
    carries the same 122,848 chromosomes and `frequency` sums to 1.0 per gene.

    NULL ALLELES ARE DROPPED and the rest renormalised. jMorp assigns every chromosome
    an allele, and at DRB4 it uses the null DRB4*03:01N — 0.5430 of all chromosomes —
    to mean "this haplotype carries no functional DRB4". That is what hla.typing
    records as `not_typed` and leaves out of its denominator, so keeping the nulls
    compares two different quantities. It is the entire DRB4 discrepancy: max|delta|
    against our controls is exactly 0.5430, and dropping them moves DRB4 from r = 0.625
    to r = 0.991, with the carriage rates then agreeing independently (ours 0.4018,
    jMorp non-null 0.4565).

    DRB3 does NOT reconcile this way and is not made to. It has no null alleles at all
    yet still spans all 122,848 chromosomes, and its DRB3*01:01P is 0.674 against our
    0.180. jMorp does not document the rule it applies there; see
    docs/OPEN_QUESTIONS.md.
    """
    d = pd.read_csv(truth_path)
    for col in ('gene', 'allele', 'count'):
        if col not in d.columns:
            raise SystemExit(f'ABORT: {truth_path} has no {col!r} column. '
                             f'Found: {list(d.columns)}')
    out, dropped = {}, {}
    for g in genes:
        sub = d[d.gene == g]
        if not len(sub):
            # A gene the panel does not carry is recorded as having no reference,
            # never silently compared against nothing.
            out[g] = ({}, 0)
            dropped[g] = []
            continue
        counts, n, nulls = {}, 0, []
        for _, r in sub.iterrows():
            a = str(r['allele'])
            name = a.split('*', 1)[1] if '*' in a else a
            if name.rstrip('P').endswith('N'):
                nulls.append(a)
                continue
            k = to_pgroup(f'{g}*{name}', pmap) if not name.endswith('P') else f'{g}*{name}'
            c = int(r['count'])
            counts[k] = counts.get(k, 0) + c
            n += c
        out[g] = (counts, n)
        dropped[g] = nulls
    return out, dropped


def reference_counts(truth_path, population, genes, depth, pmap=None):
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
            key = to_pgroup(f'{g}*{trim(m.group(1), depth)}', pmap)
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

    # This comparison is ancestry-sensitive: it asks whether an allele we called
    # appears in a mainland Japanese panel, so a control that is not mainland
    # Japanese fails it for a reason that has nothing to do with typing. The
    # restriction that makes it valid is applied ONCE, upstream, by
    # build_manifest.py --cohort-keep, and every table this component publishes
    # carries it. There is deliberately no second restriction here; two mechanisms
    # for one rule is how the numbers start disagreeing.
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

    if args.truth_format == 'jmorp_long' and not args.p_groups:
        raise SystemExit('ABORT: --truth-format jmorp_long needs --p-groups. jMorp names '
                         'alleles by P group and we do not; without the translation the '
                         'comparison is silently wrong, not merely coarser.')
    pmap = load_pgroups(args.p_groups) if args.p_groups else None

    o_ctrl = observed_counts(ctrl, genes, args.field_depth, pmap)
    o_case = observed_counts(case, genes, args.field_depth, pmap)
    if args.truth_format == 'jmorp_long':
        o_ref, ref_dropped = reference_counts_jmorp(args.truth, genes, pmap)
        n_dropped = sum(len(v) for v in ref_dropped.values())
        if n_dropped:
            print(f'[allele_freq_check] {n_dropped} null (N) reference allele(s) dropped '
                  f'and the rest renormalised — a null means the gene is absent from the '
                  f'haplotype, which this pipeline records as `not_typed` and leaves out '
                  f'of its denominator:', file=sys.stderr)
            for g, vals in ref_dropped.items():
                if vals:
                    print(f'    {g}: {", ".join(sorted(vals))}', file=sys.stderr)
    else:
        o_ref, ref_dropped = reference_counts(args.truth, args.population, genes,
                                              args.field_depth, pmap)
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
            # reference is itself an estimate from a finite panel, and treating it
            # as known would make every difference look significant. This script
            # runs against panels three orders of magnitude apart in size, so the
            # panel's own uncertainty cannot be assumed away for either.
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

    # ── per sample, so the confound has something to be plotted against ──────
    # An allele the reference never carries is not necessarily wrong: a panel
    # samples a finite number of chromosomes and cannot contain an allele rarer
    # than about one over that number. How much that matters is a property of the
    # PANEL, which is why both are run. The SHARE of such calls is still the only
    # per-sample quantity here that tracks typing quality rather than completeness.
    confirmed = {g: set(o_ref[g][0]) for g in genes}
    srows = []
    for _, r in calls.iterrows():
        n = b = 0
        for g in genes:
            if r.get(f'{g}_state', '') not in USABLE:
                continue
            for a in (r.get(f'{g}_1', ''), r.get(f'{g}_2', '')):
                if not a:
                    continue
                n += 1
                # Must go through the SAME translation the reference keys did,
                # or every P-group-named reference allele reads as unconfirmed.
                b += to_pgroup(trim(a, args.field_depth), pmap) not in confirmed[g]
        srows.append({'sample_id': r['sample_id'], 'group': r['__group'],
                      'n_chr': n, 'n_unconfirmed': b,
                      'frac_unconfirmed': round(b / n, 5) if n else ''})
    sample_tab = pd.DataFrame(srows, columns=['sample_id', 'group', 'n_chr',
                                              'n_unconfirmed', 'frac_unconfirmed'])
    if 'platform' in man.columns:
        sample_tab = sample_tab.merge(
            man[[c for c in ('sample_id', 'platform', 'observed_depth') if c in man.columns]],
            on='sample_id', how='left')
    sample_tab.to_csv(args.out_sample, sep='\t', index=False)

    cols = ['gene', 'allele', 'count_ctrl', 'n_chr_ctrl', 'freq_ctrl', 'ci_lo_ctrl',
            'ci_hi_ctrl', 'count_case', 'n_chr_case', 'freq_case', 'count_ref',
            'n_chr_ref', 'freq_ref', 'diff_ctrl_ref', 'fisher_p_ctrl_ref',
            'in_ref_only', 'in_obs_only']
    pd.DataFrame(rows, columns=cols).to_csv(args.out, sep='\t', index=False)
    scols = ['gene', 'n_alleles', 'n_chr_ctrl', 'n_chr_ref', 'pearson_r',
             'spearman_rho', 'max_abs_diff', 'allele_at_max', 'n_ref_only',
             'n_obs_only']
    pd.DataFrame(summary, columns=scols).to_csv(args.out_summary, sep='\t', index=False)

    # The reference's own denominator, reported rather than assumed: jMorp's is
    # uniform across loci, 1KG's is not (its DQB1 is typed on 186 of 210 chromosomes).
    ref_name = ('jMorp 61KJPN-HLA' if args.truth_format == 'jmorp_long'
                else f'1000 Genomes {args.population}')
    ref_n = max(n for _, n in o_ref.values()) // 2
    print(f'[allele_freq_check] {len(ctrl)} control(s) vs {ref_name} '
          f'(up to {ref_n:,} individuals) over {len(genes)} locus/loci')
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
    # ---- the crosswalk the association reads -------------------------------
    # Covers EVERY gene this cohort carries calls for, not just the reference loci:
    # a hit at DMA has to be lookup-able too, it just comes back in_reference = 0.
    all_genes = sorted({c[:-2] for c in calls.columns if c.endswith('_1')})
    ref_freq = {}
    for g in genes:
        cr, nr = o_ref[g]
        if nr:
            ref_freq.update({k: v / nr for k, v in cr.items()})
    pg_rows = []
    for g in all_genes:
        c1, c2, cs = f'{g}_1', f'{g}_2', f'{g}_state'
        seen = set()
        for _, r in calls.iterrows():
            if r.get(cs, '') not in USABLE:
                continue
            for a in (r.get(c1, ''), r.get(c2, '')):
                if a:
                    seen.add(trim(a, args.field_depth))
        for k in sorted(seen):
            pg = to_pgroup(k, pmap)
            pg_rows.append({'gene': g, 'allele': k, 'p_group': pg,
                            'in_reference': int(pg in ref_freq),
                            'freq_reference': round(ref_freq[pg], 6) if pg in ref_freq else ''})
    pg_tab = pd.DataFrame(pg_rows, columns=['gene', 'allele', 'p_group',
                                            'in_reference', 'freq_reference'])
    pg_tab.to_csv(args.out_pgroup_map, sep='\t', index=False)

    print(f'    {len(rows)} (gene, allele) row(s) -> {args.out}')
    n_abs = int((pg_tab.in_reference == 0).sum())
    print(f'    {len(pg_tab):,} allele(s) over {len(all_genes)} gene(s) -> '
          f'{args.out_pgroup_map}; {n_abs:,} ({n_abs / len(pg_tab):.1%}) of the DISTINCT '
          f'alleles are absent from the reference — read that beside the '
          f'chromosome-weighted figure below, which is the one that matters: most absent '
          f'alleles are carried by one or two chromosomes each. in_reference = 0 is a '
          f'per-hit artefact flag for the association, not an error.')
    tot_n, tot_b = int(sample_tab.n_chr.sum()), int(sample_tab.n_unconfirmed.sum())
    print(f'    {len(sample_tab):,} sample row(s), {tot_b:,}/{tot_n:,} '
          f'({tot_b / tot_n:.1%}) chromosomes on an allele the reference does not carry '
          f'-> {args.out_sample}')
    print('    This compares FREQUENCIES. It does not measure per-sample accuracy.')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in allele_freq_check: {e}', file=sys.stderr)
        sys.exit(1)
