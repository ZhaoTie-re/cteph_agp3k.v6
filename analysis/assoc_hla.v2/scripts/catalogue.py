#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : ONE document per figure family (results/figures/0N.<family>/README.md).
#
#           A family's prose -- the question, what each panel shows, how to read
#           it, what it cannot establish -- is written here once. Every number in
#           it is computed from the same tables the figures read (figures.Data)
#           or from the members' <name>.stats.json, so the document cannot
#           disagree with the PNGs beside it.
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_doc                      # noqa: E402
import vocab as V                      # noqa: E402
from figures import (COMPARISON, FAMILIES, PM, Data, add_data_args, majority,   # noqa: E402
                     nlp)

METHODS_REF = '../../../docs/METHODS.md'     # results/figures/<family>/README.md -> docs/
SUPPLEMENTARY = {'thresholds', 'models', 'qq'}
COMMON_LIMITS = [
    'Replication: the three cohorts are nested and share almost all cases; n/3 is a sample-selection '
    'sensitivity analysis, not independent replication.',
    'Platform: all cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; a typing difference '
    'between platforms would look like an association, and no covariate separates the two. Within cases the '
    'platform mix also differs by DVT status.',
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--family', required=True, choices=FAMILIES)
    add_data_args(p)
    p.add_argument('--stats', nargs='*', default=[])
    p.add_argument('--figure-dir', required=True)
    p.add_argument('--out', default='README.md')
    return p.parse_args()


def members(paths):
    out = {}
    for p in paths:
        d = json.load(open(p))
        out[d['peak']] = d['values']
    return dict(sorted(out.items()))


def fmt3(xs, f='{:.1e}'):
    return ' / '.join('-' if pd.isna(x) else f.format(x) for x in xs)


def signal_text(D, cs=None):
    s = D.signals if cs is None else D.signals[D.signals.case_set == cs]
    return '; '.join(f"{V.marker_label(r.lead, r[3])} ({V.CASE_SET_SHORT[r.case_set]}, {r[3]}: GEC {r.gec}, "
                     f"max-T {r.maxt}, Bonferroni {r.bonf})" for r in s.itertuples()) or 'none'


# ═══════════════════════════════════════════════════════════════ per family
def doc_landscape(D, mem):
    rows = []
    for cs in V.CASE_SETS:
        m = D.summ[D.summ.case_set == cs]
        for cls in V.CLASSES:
            tests = set(D.tests(cs, cls))
            e = m[(m['class'] == cls) & m.marker.isin(tests)]
            rec = {'case set': V.CASE_SET_SHORT[cs], 'class': cls, 'tests drawn': len(tests)}
            for t in V.THRESHOLDS:
                rec[f'{t}: 3/3 / 2/3 / 1/3'] = ' / '.join(str(int((e[f'{t}_n'] == k).sum())) for k in (3, 2, 1))
            rows.append(rec)
    n_dq = int(((D.summ.case_set == 'dvt_no') & (D.summ['class'] == 'residue') & (D.summ.gec_n > 0)).sum())
    interp = (f'Random model. Independent signals (02.signals, 03.conditioning): {signal_text(D)}. In the bar '
              'encoding, a test whose whole bar lies above the solid line is significant in all three cohorts; a '
              'bar that crosses it with its top only is significant in one. Residues on one haplotype rise '
              f'together, so a run of coloured bars can be one signal: in DVT no, {n_dq} residue markers are '
              'significant in some cohort, and the conditional analysis (03.conditioning) asks how many signals '
              'they are.')
    panels = [('bars.<case_set>.png a-b', 'one bar per test, alleles (a) and residues (b)',
               'For each test (a residue mirror pair is drawn once) the bar runs from −log10 of its largest P to '
               '−log10 of its smallest P over the three cohorts, random model: its top is the strongest cohort and '
               'its bottom the weakest. Colour and line width give the number of cohorts in which the test is '
               'significant at that cohort\'s LD-adjusted (GEC) threshold with the same direction (3/3 darkest, '
               'then 2/3, 1/3 "Single cohort"; grey not significant). Horizontal lines: GEC (solid), exact max-T '
               '(dashed) and conventional Bonferroni (dotted), each the median of the three cohorts\' thresholds '
               '(they differ by < 1 %). The violet triangle marks where the genome-wide MHC lead lies among the '
               'genes. Signal leads are labelled with their n/3.'),
              ('manhattan.<case_set>.png a-f', 'every test in each cohort',
               'Rows narrow, intermediate, full mainland; columns alleles, residues. One point per test, red if '
               'significant at that cohort\'s own GEC threshold; the three lines are the cohort\'s own thresholds. '
               'The same data as the bars, cohort by cohort.')]
    reading = ['Start with bars.all.png: find the bars above the solid line and their colour (n/3).',
               'Check the dashed and dotted lines: does the call survive the exact and the conventional threshold?',
               'Open manhattan.<case_set>.png of the same case set to see which cohort carries a 1/3 call.',
               'Compare case sets (bars.dvt_yes / dvt_no / dvt_case_only) with 02.signals for the effect sizes.']
    limits = ['Independence of coloured bars: neighbouring residues and alleles share haplotypes; see '
              '03.conditioning.', 'Physical position: the x axis orders genes and markers; it is not a coordinate.']
    return ('Where in the MHC, and in how many cohorts',
            'Along the typed HLA genes, which allele and residue tests are associated in each case set, and in how '
            'many of the three nested cohorts does each hold?', panels, interp,
            [('Tests by cross-cohort consistency (random model), per case set and class', pd.DataFrame(rows))],
            reading, limits)


def doc_signals(D, mem):
    s = D.signals
    rows = []
    for r in s.itertuples():
        cls = r[3]
        for comp in V.CASE_SETS:
            e = D.primary(comp, cls)
            e = e[e.marker == r.lead].set_index('cohort').reindex(list(V.COHORTS))
            rec = {'signal': V.marker_label(r.lead, cls), 'selected in': V.CASE_SET_SHORT[r.case_set],
                   'comparison': COMPARISON[comp], 'OR N / I / F': fmt3(e['or'], '{:.2f}'),
                   'P N / I / F': fmt3(e['p'])}
            for t in V.THRESHOLDS:
                rec[t] = f"{majority(e[f'sig_{t}'].fillna(False).astype(bool).tolist(), e['or'].tolist())}/3"
            rows.append(rec)
    sig = s[['case_set', 'class', 'signal', 'lead', 'gec', 'maxt', 'bonf', 'n_members', 'carrier_alleles']].copy()
    sig['case_set'] = sig.case_set.map(V.CASE_SET_SHORT)
    interp = (f'{len(s)} independent signal(s) in the random model: {signal_text(D)}. Each is shown in every '
              'comparison, so a signal selected in one case set can be read in the others. An effect of the same '
              'size and direction in DVT yes and DVT no, with a case-only OR near 1, is a shared effect that the '
              'smaller DVT sets lack the power to call; a case-only OR away from 1, significant, would be the '
              'evidence of a difference. The DVT-yes and DVT-no comparisons share their controls, so their '
              'side-by-side estimates are not a test.')
    panels = [('signals.png a', 'odds ratio in every comparison and cohort',
               'Rows grouped by signal (bold: the lead, and the case set it was selected in); within a group, one '
               'row per comparison: all cases, DVT yes and DVT no against all controls, and DVT yes against DVT '
               'no within cases. Three estimates per row, one per cohort (dark narrow → light full), OR per copy '
               'with 95 % CI, random model. Filled: significant at that comparison\'s LD-adjusted threshold.'),
              ('signals.png b', 'significance at each threshold',
               'For the same rows: significance in each cohort (N, I, F) at the LD-adjusted (GEC), exact (max-T) '
               'and conventional Bonferroni thresholds of that comparison, cohort and class (filled = '
               'significant), and the number of cohorts with the same direction (n/3).')]
    reading = ['Read the first row of each group: the comparison the signal was selected in.',
               'Read down the group: is the OR the same in DVT yes and DVT no, and is the case-only OR near 1?',
               'Read panel b across: does the call hold at max-T and at conventional Bonferroni?']
    limits = ['Fine-mapping: the lead stands for the markers it explains (03.conditioning, members).',
              'Random-model ORs are SAIGE one-step estimates; read them for direction and rough size.']
    return ('The independent signals', 'For each independent signal, how large is the effect in every case-set '
            'comparison and cohort, and at which thresholds is it significant?', panels, interp,
            [('Independent signals (random model)', sig), ('Every signal in every comparison', pd.DataFrame(rows))],
            reading, limits)


def doc_conditioning(D, mem):
    st = D.steps[D.steps.model == PM].copy()
    st['case_set'] = st.case_set.map(V.CASE_SET_SHORT)
    st = st[['case_set', 'class', 'step', 'conditioned_on', 'lead', 'entered', 'reason', 'gec', 'n_sig_gec_any',
             'p_narrow', 'p_intermediate', 'p_full']]
    car = []
    for r in D.signals[D.signals['class'] == 'residue'].itertuples():
        gene = V.parse_marker(r.lead, 'residue')[0]
        carriers = {x.strip() for x in str(r.carrier_alleles).split(',') if x.strip()}
        al = D.primary(r.case_set, 'allele')
        for m in sorted({x for x in al.marker if x.startswith(f'{gene}*')}, key=lambda n: V.sort_key(n, 'allele')):
            e = al[al.marker == m].set_index('cohort').reindex(list(V.COHORTS))
            fc = e.loc['full_mainland', 'freq_control']
            car.append({'signal': V.marker_label(r.lead, 'residue'), 'case set': V.CASE_SET_SHORT[r.case_set],
                        'allele': m, 'carries the residue': m in carriers,
                        'OR N / I / F': fmt3(e['or'], '{:.2f}'), 'P N / I / F': fmt3(e['p']),
                        'control frequency (full)': round(float(fc), 4) if pd.notna(fc) else np.nan})
    left = {}
    for (cs, cls), g in D.steps[D.steps.model == PM].groupby(['case_set', 'class']):
        if g.entered.any():
            left[(cs, cls)] = int(g.loc[g.step == g.step.max(), 'n_sig_gec_any'].iloc[0])
    interp = (f'Random model. Signals: {signal_text(D)}. After conditioning on the lead(s), the number of tests '
              'still significant at the LD-adjusted threshold in some cohort is ' +
              ', '.join(f'{V.CASE_SET_SHORT[c]} {k}: {n}' for (c, k), n in left.items()) +
              '. A walk that enters one signal and leaves nothing explains all its significant markers with one '
              'lead. The carrier panels show which tested alleles carry each residue signal: a residue-level '
              'explanation predicts that the carriers share the residue\'s direction and the other alleles sit '
              'near OR 1.')
    panels = [('conditioning.png', 'each walk before and after conditioning',
               'One row per (case set, class) walk that entered a signal, random model. Left: the marginal scan '
               '(step 0); right: the scan conditional on every entered lead (last step). Bars as in 01.landscape '
               '(largest to smallest P over the cohorts, coloured by n/3 at GEC). The triangle on the right marks '
               'the conditioned lead.'),
              ('carriers.png', 'the alleles that carry each residue signal',
               'One panel per residue signal: every tested allele of the gene, OR per copy with 95 % CI in each '
               'cohort (random model, the signal\'s own comparison), carriers of the residue first. Carriers are '
               'read from hla.typing.v2: allele_match_map (the reference allele a 2-field allele was matched to) '
               'and residue_reference (that allele\'s residue at the position).')]
    reading = ['conditioning.png: compare left and right in each row; the coloured bars should turn grey.',
               'carriers.png: do the carriers share the residue\'s direction, and do the non-carriers sit near 1?']
    limits = ['Statistical, not functional, separation: markers in near-perfect LD with the lead cannot be told '
              'apart from it.', 'Untested alleles: only alleles that passed the typing QC (MAF ≥ 0.5 %) are shown.']
    return ('Conditional analyses: independent signals and their carriers',
            'Does one marker explain each set of significant markers, and which alleles carry each residue '
            'signal?', panels, interp,
            [('Every step of every walk (random model)', st),
             ('Alleles of each residue signal\'s gene', pd.DataFrame(car))], reading, limits)


def doc_gwas_link(D, mem):
    if D.gwas.empty:
        return ('The GWAS MHC signal', 'No GWAS MHC lead passed the three-cohort rule.', [], None, [], [], [])
    ldid = D.gwas.lead_id.iloc[0]
    lg = D.lg[(D.lg.case_set == 'all') & (D.lg.model == PM) & (D.lg.lead == ldid)
              & ~D.lg.why.str.contains('jointly')]
    w = lg.pivot_table(index=['conditioned_on', 'why'], columns='cohort', values='lead_p_given', aggfunc='first')
    w = w.reindex(columns=list(V.COHORTS))
    r2 = lg[lg.cohort == V.COHORTS[0]].groupby('conditioned_on').r.first() ** 2
    tab = pd.DataFrame({'conditioning marker': [i[0] for i in w.index], 'why': [i[1] for i in w.index],
                        'r2 with lead (narrow)': [round(float(r2.get(i[0], np.nan)), 3) for i in w.index],
                        'lead P given it N / I / F': [fmt3(w.loc[i]) for i in w.index]})
    tab = tab.sort_values('r2 with lead (narrow)', ascending=False)
    marg = lg.groupby('cohort').lead_p_marginal.first().reindex(list(V.COHORTS))
    vals = [(k, v) for vs in mem.values() for k, v in vs]
    interp = (f'All cases, random model. The GWAS lead {ldid} has marginal P {fmt3(marg)} here (mean-imputed; the '
              'genome-wide scan, with best_guess imputation, reports a slightly smaller P, METHODS §7.1). '
              + ' '.join(f'{k}: {v}.' for k, v in vals) +
              ' A marker that explains the lead pulls its P towards 1; one in LD that does not leaves it near the '
              'marginal line.')
    panels = [('gwas_link.png a', 'the GWAS lead per cohort',
               'Marginal −log10 P (black), given the first residue signal\'s lead (red), and given the HLA marker '
               'in highest LD with it (open square); dotted: the GWAS suggestive line 1e-5.'),
              ('gwas_link.png b', 'the HLA signal lead per cohort',
               'Marginal (filled) and given the GWAS lead (open); tick: the cohort\'s residue GEC threshold.'),
              ('gwas_link.png c-e', 'every single-marker conditioning, per cohort',
               'x: r² of the conditioning marker with the GWAS lead (ALT dosage); y: the lead\'s −log10 P given that '
               'marker. Markers tried: the 10 in highest |r|, the 10 with the smallest marginal P and the signal '
               'leads. Dashed: the lead\'s marginal −log10 P.')]
    reading = ['Panel a: does conditioning on the HLA signal pull the lead down in every cohort, and the LD tag not?',
               'Panel b: does the HLA signal lose its significance given the lead (the reverse)?',
               'Panels c-e: is the marker that explains the lead the one with the highest r²?']
    limits = ['Which variant is causal: two markers that explain each other are one signal statistically.',
              'The lead is a SNP between HLA-DRB1 and HLA-DQA1; its relation to DQB1 runs through haplotypes.']
    return ('The genome-wide MHC signal and the HLA signal',
            'Is the genome-wide scan\'s MHC lead the same signal as an HLA residue or allele?', panels, interp,
            [('The GWAS lead given each single HLA marker (all cases, random model)', tab)], reading, limits)


def doc_dvt(D, mem):
    rows = []
    for m, cls in sorted({(r.lead, r[3]) for r in D.signals.itertuples()}):
        for coh in V.COHORTS:
            rec = {'marker': V.marker_label(m, cls), 'cohort': V.COHORT_SHORT[coh]}
            for cs in V.CASE_SETS:
                e = D.primary(cs, cls, coh)
                e = e[e.marker == m]
                rec[f'OR {COMPARISON[cs]}'] = round(float(e['or'].iloc[0]), 3) if len(e) else np.nan
                rec[f'P {COMPARISON[cs]}'] = float(e['p'].iloc[0]) if len(e) else np.nan
            rows.append(rec)
    rr = [float(str(v).split(', ')[1]) for vs in mem.values() for k, v in vs if 'r(log OR)' in k]
    interp = ('Random model, one point per test. r(log OR) between the two sets is '
              f'{min(rr):.2f}-{max(rr):.2f} over the six panels. Points far from the identity line differ between '
              'cases with and without DVT, but the two sets share their controls, so the scatter is not a test; '
              'the case-only test (rings) is. The DVT sets are small (140-159 cases), so a marker significant in '
              'one set and not the other is expected by chance for some markers.') if rr else None
    panels = [('dvt.png a-f', 'log OR with DVT against without DVT',
               'Rows alleles, residues; columns narrow, intermediate, full mainland. x: log OR of cases without DVT '
               'against all controls; y: cases with DVT against the same controls; random model, one point per '
               'test. Blue / red: significant at the LD-adjusted threshold in DVT no / DVT yes; ring: significant '
               'in the case-only test. Grey line: identity. Signal leads are labelled.')]
    reading = ['Points on the identity line: the same effect with and without DVT.',
               'A coloured point off the line: significant in one set with a different estimate in the other.',
               'A ring: the case-only test supports a difference.']
    limits = ['DVT status is recorded only for cases; about 30 % of cases have unknown status and are in neither set.',
              'A difference read off the scatter alone: the two axes share controls; use the case-only test.']
    return ('Cases with and without DVT', 'Does any HLA allele or residue associate differently with CTEPH in cases '
            'with a history of DVT than in cases without?', panels, interp,
            [('Signal leads in every comparison (random model)', pd.DataFrame(rows))], reading, limits)


def doc_thresholds(D, mem):
    t = D.thr.copy()
    t['case_set'] = t.case_set.map(V.CASE_SET_SHORT)
    t['cohort'] = t.cohort.map(V.COHORT_SHORT)
    tab = t[['case_set', 'cohort', 'class', 'markers', 'mirror_pairs', 'tests', 'me_gec', 'me_maxt', 'p_gec',
             'p_maxt', 'p_bonf', 'fwer_gec', 'fwer_bonf']].round({'me_gec': 1, 'me_maxt': 1, 'fwer_gec': 3,
                                                                  'fwer_bonf': 3})
    al, rs = D.thr[D.thr['class'] == 'allele'], D.thr[D.thr['class'] == 'residue']
    interp = (f'Over the 24 (case set, cohort, class) families the GEC threshold\'s family-wise error under max-T is '
              f'{al.fwer_gec.min():.3f}-{al.fwer_gec.max():.3f} for alleles and '
              f'{rs.fwer_gec.min():.3f}-{rs.fwer_gec.max():.3f} for residues; GEC v0.2 estimates Me within blocks '
              'of at most 100 markers and ignores LD between blocks. Conventional Bonferroni counts every test as '
              f'independent: {al.fwer_bonf.min():.3f}-{al.fwer_bonf.max():.3f} for alleles and '
              f'{rs.fwer_bonf.min():.3f}-{rs.fwer_bonf.max():.3f} for residues, whose mirror pairs and haplotypes '
              'make them strongly correlated.')
    panels = [('thresholds.png a', 'alleles: tests and effective numbers',
               'Per case set and cohort: the number of tests (conventional Bonferroni divides by it), the Me '
               'implied by the exact max-T threshold (0.05 / threshold), and GEC v0.2\'s Me (the LD-adjusted '
               'threshold divides by it).'),
              ('thresholds.png b', 'residues: the same', 'A residue mirror pair counts as one test.'),
              ('thresholds.png c', 'family-wise error rate actually achieved',
               'For each of the 24 families, the probability under the max-T simulation (400,000 multivariate-normal '
               'draws on the covariate-adjusted genotype correlation) that any test passes the GEC (filled) or the '
               'conventional Bonferroni (open) threshold; dashed: the nominal 0.05.')]
    reading = ['Compare the red dot (GEC) with the grey diamond (exact): the closer, the better calibrated.',
               'Panel c: points on the dashed line hold the error rate exactly; below it, the threshold is stricter '
               'than needed.']
    limits = ['The max-T simulation assumes the test statistics are jointly normal with the genotype correlation; '
              'it is a calibration, not a permutation of phenotypes.']
    return ('Supplementary: the three thresholds', 'How strict is each threshold, and how close does each come to '
            'a family-wise error rate of 0.05?', panels, interp, [('Every family', tab)], reading, limits)


def doc_models(D, mem):
    rows = []
    for cs in V.CASE_SETS:
        for coh in V.COHORTS:
            for cls in V.CLASSES:
                d = D.assoc[(D.assoc.case_set == cs) & (D.assoc.cohort == coh) & (D.assoc['class'] == cls)
                            & (D.assoc.test_id == D.assoc.marker)]
                p = d.pivot_table(index='marker', columns='model', values='p', aggfunc='first').dropna()
                b = d.pivot_table(index='marker', columns='model', values='beta', aggfunc='first').reindex(p.index)
                g = d.pivot_table(index='marker', columns='model', values='sig_gec',
                                  aggfunc='first').fillna(False).astype(bool)
                rows.append({'case set': V.CASE_SET_SHORT[cs], 'cohort': V.COHORT_SHORT[coh], 'class': cls,
                             'tests': len(p),
                             'r(-log10 P)': round(float(np.corrcoef(nlp(p['fixed']), nlp(p['random']))[0, 1]), 4),
                             'r(beta)': round(float(np.corrcoef(b['fixed'], b['random'])[0, 1]), 4),
                             'GEC fixed / random / both': f"{int(g['fixed'].sum())} / {int(g['random'].sum())} / "
                                                          f"{int((g['fixed'] & g['random']).sum())}"})
    tab = pd.DataFrame(rows)
    vals = dict((k, v) for vs in mem.values() for k, v in vs)
    interp = ('The two models differ in samples (the fixed model drops relatives), missing calls (complete case '
              'against mean imputation) and estimator (Wald against score), so points scatter around the identity; '
              f"r(−log10 P) is {tab['r(-log10 P)'].min():.3f}-{tab['r(-log10 P)'].max():.3f} per case set, cohort "
              'and class. Tests significant at the LD-adjusted threshold in one model only: '
              f"{vals.get('tests significant at GEC in one model only (cohort, case set, marker)', 'none')}. "
              'The random model is primary; calls made by the fixed model alone are not counted.')
    panels = [('models.png a-c', '−log10 P, fixed (x) against random (y)', 'One panel per cohort; one point per '
               'test, all four case sets (colour). Grey: identity.'),
              ('models.png d-f', 'β, fixed against random', 'Fixed: log OR from plink2. Random: SAIGE\'s one-step '
               'estimate Tstat / var.')]
    reading = ['Points on the identity: the models agree.', 'The table and the interpretation list the calls '
               'that differ.']
    limits = ['Which model is right for a marker: they answer the same question under different assumptions.']
    return ('Supplementary: fixed against random model', 'Does accounting for relatedness and the polygenic '
            'background (random model) change any test compared with the fixed model?', panels, interp,
            [('Per case set, cohort and class', tab)], reading, limits)


def doc_qq(D, mem):
    vals = [(k, v) for vs in mem.values() for k, v in vs]
    lam = [float(str(v).split(', ')[1]) for _k, v in vals]
    co = [float(str(v).split(', ')[1]) for k, v in vals if k.startswith(V.CASE_SET_SHORT['dvt_case_only'])]
    interp = (f'λGC spans {min(lam):.2f}-{max(lam):.2f} over the 24 scans (random model), '
              f'{min(co):.2f}-{max(co):.2f} in the case-only test. All markers come from one region of chromosome '
              '6 and are strongly correlated, so a few real or artefactual differences between the groups lift the '
              'whole curve; in the genome-wide fixed scan the MHC SNPs themselves have λGC 1.89 against 1.13 '
              'elsewhere (docs/STUDY_NOTES.md §8). λGC here describes the region, not the calibration of the test, '
              'and is never applied.')
    panels = [('qq.png a-h', 'observed against expected −log10 P',
               'Rows case sets, columns classes; one curve per cohort (random model, one point per test); grey: '
               'pointwise 95 % band of the uniform null, Beta(i, n − i + 1); corner: λGC per cohort.')]
    reading = ['Compare the three curves in a panel: the cohorts should behave alike.',
               'Do not read λGC as genome-wide inflation.']
    limits = ['Whether the departure is biology or artefact: see the platform limit below.']
    return ('Supplementary: P values against the null', 'How far does each scan depart from the uniform null?',
            panels, interp, [('λGC per scan', pd.DataFrame(vals, columns=['scan', 'tests, lambda_GC']))],
            reading, limits)


def main():
    a = parse_args()
    D = Data(a)
    mem = members(a.stats)
    title, question, panels, interp, tables, reading, limits = globals()[f'doc_{a.family}'](D, mem)
    role = 'Supplementary figure' if a.family in SUPPLEMENTARY else 'Main figure'
    subject = f"**{role}.** **Covers:** " + ', '.join(f'`{k}.png`' for k in mem) + f' in `{a.figure_dir}/`'
    numbers = [(f'{k}: {lab}', v) for k, vals in mem.items() for lab, v in vals]
    figure_doc.write_doc(None, title=title, question=question, panels=panels, interpretation=interp,
                         numbers=numbers, tables=tables, reading=reading, limits=limits + COMMON_LIMITS,
                         methods_ref=METHODS_REF, out_path=a.out, subject=subject)
    print(f'[catalogue] {a.family}: {len(mem)} member(s) -> {a.out}')


if __name__ == '__main__':
    main()
