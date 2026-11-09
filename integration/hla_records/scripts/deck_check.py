#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Agreement of the workbook's tables with the two decks' data. Every
#           quantity the workbook shares with a slide is compared with the value
#           the deck's build read: the typing deck's stats.json; the association
#           deck's stats.json, data/results_all.json and data/pos57_pgroups.tsv.
#           The map is explicit, sheet by sheet; counts must be equal, floats
#           agree to a relative 1e-6 (the tables hold 7 significant digits; the
#           genome-wide lambda 1e-5, its P values are read at another precision),
#           text must be equal. A mismatch is never resolved by copying the deck:
#           it fails the run and is listed for tracing.
# Inputs  : --typing --assoc (the table folders) --typing-stats --assoc-stats
#           --assoc-results --assoc-pgroups
# Outputs : --out (deck_check.json: one check per sheet, every mismatch listed)
# Component: hla_records
# ---------------------------------------------------------------------------
import argparse
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import pandas as pd

COH = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
SHORT = {'narrow_mainland': 'narrow', 'intermediate_mainland': 'intermediate', 'full_mainland': 'full'}


class Cmp:
    def __init__(self):
        self.n = defaultdict(int)
        self.bad = defaultdict(list)

    def __call__(self, sheet, item, ours, deck, rel=1e-6):
        self.n[sheet] += 1
        ok = same(ours, deck, rel)
        if not ok:
            self.bad[sheet].append({'item': item, 'workbook': str(ours), 'deck': str(deck)})


def num(x):
    if isinstance(x, bool):
        return float(x)
    if isinstance(x, (int, float)):
        return float(x)
    try:
        return float(str(x).replace(',', ''))
    except ValueError:
        return None


def same(a, b, rel):
    if a is None or b is None or (isinstance(a, float) and math.isnan(a)) or (isinstance(b, float) and math.isnan(b)):
        return (a is None or (isinstance(a, float) and math.isnan(a))) and \
               (b is None or (isinstance(b, float) and math.isnan(b)))
    x, y = num(a), num(b)
    if x is not None and y is not None:
        return abs(x - y) <= rel * max(abs(x), abs(y)) + 1e-12
    return str(a) == str(b)


class Row(dict):
    """A table row with attribute access (keyword-named columns such as 'or' included)."""
    __getattr__ = dict.__getitem__


def recs(df):
    return [Row(r) for r in df.to_dict('records')]


def path(d, p):
    for k in p.split('.'):
        d = d[k]
    return d


def tsv(folder, name, dtype=None):
    return pd.read_csv(Path(folder) / name, sep='\t',
                       dtype=dtype or {'position': str, 'residue': str, 'marker': str, 'unit': str})


# ── typing deck ─────────────────────────────────────────────────────────────
def typing(c, T, S):
    for r in recs(tsv(T, 't1_run.tsv', dtype=str)):
        v = path(S, r.deck_key)
        c('T1', r.item, r.value, len(v) if isinstance(v, list) else v)
    pl = tsv(T, 't1_platforms.tsv')
    allg = pl[pl.platform == 'All genomes'].iloc[0]
    c('T1', 'outliers', allg.outliers, S['samples']['outliers'])
    for k, col in (('read_pairs<p1', 'out_read_pairs'), ('no_2field>=3', 'out_no2field'),
                   ('failed_obligate', 'out_failed_obligate')):
        c('T1', f'outliers {k}', allg[col], S['samples']['reasons'][k])
    c('T1', 'calls', allg.calls, S['calls_reported'])
    # T2
    for r in recs(tsv(T, 't2_drb.tsv')):
        d, m = S['drb'][r.gene], S['drb_matrix'][r.gene]
        for k, v in (('agree', d['agree']), ('two_where_one', d['two_where_one']),
                     ('present_where_none', d['present_where_none']), ('no_flag', d['flags']['']),
                     ('reconciled', d['flags'].get('drb_reconciled', 0)), ('unexpected', d['flags'].get('drb_unexpected', 0)),
                     ('missing', d['flags'].get('drb_missing', 0)), ('drb1_unknown', d['flags'].get('drb1_unknown', 0)),
                     ('with_drb1', m['with_drb1']), ('with_result', m['with_result']), ('exported_match', m['match']),
                     ('exported_blank', m['blank'])):
            c('T2', f'{r.gene} {k}', getattr(r, k), v)
    for r in recs(tsv(T, 't2_two_allele.tsv')):
        d = S['paralogue_two'][str(r.expected_copies)]
        c('T2', f'two-allele {r.expected_copies}', r.two_allele_calls, d['n'])
        c('T2', f'minor<15% {r.expected_copies}', r.minor_below_15pct, d['below'])
    # T3
    for r in recs(tsv(T, 't3_genes.tsv')):
        g, t = S['genes'][r.gene], S['tiers'][r.gene]
        for k in ('heterozygous', 'homozygous', 'present', 'absent', 'ambiguous'):
            c('T3', f'{r.gene} {k}', getattr(r, k), g[f'state_{k}'])
        for k in ('basis_agreed', 'basis_prior', 'basis_pgroup', 'weak_prior', 'median_min_reads', 'low_support',
                  'het_assessed', 'imbalanced', 'conf_high', 'conf_medium', 'conf_low', 'conf_no_call'):
            c('T3', f'{r.gene} {k}', getattr(r, k), g[k])
        c('T3', f'{r.gene} share_high', r.share_high, g['share_high'], rel=1e-3)   # the deck keeps 4 digits
        c('T3', f'{r.gene} calls', r.calls, g['samples'])
        c('T3', f'{r.gene} tier', r.tier, t['tier'])
        c('T3', f'{r.gene} criteria', '' if pd.isna(r.criteria_failed) else r.criteria_failed, t['criteria_failed'] or '')
        c('T3', f'{r.gene} imgt', r.imgt_alleles, S['imgt_alleles'][r.gene])
        c('T3', f'{r.gene} genotype share', r.genotype_share, S['completeness']['by_gene'][r.gene])
    for k in ('agreed', 'prior', 'prior_decisive', 'pgroup', 'none'):
        c('T3', f'basis {k}', json.loads((Path(T) / 'typing_facts.json').read_text())['totals'][k], S['basis'][k])
    # T4
    for r in recs(tsv(T, 't4_panels.tsv')):
        if r.panel.startswith('ToMMo'):
            d = S['tommo'][r.gene]
            for k in ('n_samples', 'n_chr_ours', 'r', 'D', 'excess_D', 'frequency_comparable'):
                c('T4', f'ToMMo {r.gene} {k}', getattr(r, k), d[k])
            if r.gene in S['coverage']:
                v = S['coverage'][r.gene]
                c('T4', f'{r.gene} in panel', r.alleles_in_panel, v['alleles']['in_panel'])
                c('T4', f'{r.gene} absent', r.alleles_absent, v['alleles']['absent'])
                c('T4', f'{r.gene} too new', r.alleles_too_new, v['alleles']['too_new'])
                c('T4', f'{r.gene} copies absent', r.copies_absent, v['copies']['absent'])
                c('T4', f'{r.gene} two-field', r.two_field, v['two_field'])
            pm = next((x for x in S['pgroup_mapping']['per_gene'] if x['gene'] == r.gene), None)
            if pm:
                for k, v in (('ours_keys', pm['ours_keys']), ('tommo_keys', pm['tommo_keys']), ('either', pm['either'])):
                    c('T4', f'{r.gene} {k}', getattr(r, k), v)
                c('T4', f'{r.gene} r_shared', r.r_shared, pm['r_shared'])
        else:
            d = S['kg'][r.gene]
            for k in ('n_chr_panel', 'r', 'D', 'excess_D'):
                c('T4', f'1KG {r.gene} {k}', getattr(r, k), d[k])
    for r in recs(tsv(T, 't4_harmonization.tsv', dtype=str)):
        v = path(S, r.deck_key)
        v = len(v) if isinstance(v, list) else v
        if r.deck_key == 'tommo_vintage':
            v = f'{str(v)[0]}.{str(v)[1:3]}.{str(v)[3]}'
        if isinstance(v, (int, float)) or str(r.value).isdigit():
            c('T4', r.item, r.value, v)
    fo = S['freq_outlier']
    o = dict(tsv(T, 't4_outlier.tsv', dtype=str)[['item', 'value']].values)
    c('T4', 'outlier allele', o['Allele'], fo['allele'])
    c('T4', 'outlier ToMMo copies', o['ToMMo, copies'], fo['count_panel'])
    c('T4', 'outlier ToMMo freq', o['ToMMo, frequency'], fo['freq_panel'], rel=1e-4)
    c('T4', 'outlier name now', o['Name in IPD-IMGT/HLA 3.64.0'], fo['name_now'])
    c('T4', 'outlier group ours', o['DRB1*15 group, our controls'], fo['group']['ours'], rel=1e-3)
    c('T4', 'outlier group ToMMo', o['DRB1*15 group, ToMMo'], fo['group']['panel'], rel=1e-3)
    for r in recs(tsv(T, 't4_homozygosity.tsv')):
        for k in ('obs_hom', 'exp_hom', 'ratio', 'dropout', 'dropout_lo', 'dropout_hi'):
            c('T4', f'{r.gene} {k}', getattr(r, k), S['hom'][r.gene][k])
    cv = tsv(T, 't4_confidence_vs_panel.tsv')
    for r in recs(cv[cv.gene == 'all ToMMo loci']):
        for k in ('n_chr', 'n_unconfirmed', 'unconfirmed_share'):
            c('T4', f'conf {r.confidence} {k}', getattr(r, k), S['conf_vs_panel'][r.confidence][k])
    tf = json.loads((Path(T) / 'typing_facts.json').read_text())
    c('T4', 'shift n', tf['shift']['n'], S['shift']['n'])
    c('T4', 'shift ours', tf['shift']['kg_on_our_side'], S['shift']['ours'])
    c('T4', 'shift tommo', tf['shift']['kg_on_tommo_side'], S['shift']['tommo'])
    # T5
    rec = tsv(T, 't5_record.tsv')
    q = S['qc_per_gene']['per_cohort']
    for r in recs(rec):
        cls = 'alleles' if r['class'] == 'allele' else 'residues'
        c('T5', f'{r.cohort} {r.stage} {cls}', r.columns, sum(q[r.cohort][r.stage][cls].values()))
    g = tsv(T, 't5_genes.tsv')
    for r in recs(g):
        if r.gene == 'All':
            continue
        stages = (('typed', 'allele_typed', 'residue_typed'), ('callrate', 'allele_callrate', 'residue_callrate'),
                  ('maf', 'allele_maf', 'residue_maf'))
        if r.gene in ('DRB3', 'DRB4', 'DRB5'):   # not tested: after 'typed' the sheet marks them not applicable
            stages = stages[:1]
        for stage, ka, kr in stages:
            for cls, k in (('alleles', ka), ('residues', kr)):
                v = q[r.cohort][stage][cls].get(r.gene, 0)
                c('T5', f'{r.cohort} {r.gene} {stage} {cls}', getattr(r, k), v)
        if r.gene in S['allele_coverage'][r.cohort]['by_gene']:
            c('T5', f'{r.cohort} {r.gene} coverage', r.allele_copies_covered, S['allele_coverage'][r.cohort]['by_gene'][r.gene])
    for co in COH:
        a = g[(g.cohort == co) & (g.gene == 'All')].iloc[0]
        c('T5', f'{co} coverage all', a.allele_copies_covered, S['allele_coverage'][co]['all'])
        c('T5', f'{co} residue coverage', a.residue_copies_covered, S['residue_coverage'][co]['all'])
        c('T5', f'{co} positions kept', a.residue_positions_kept, S['residue_coverage'][co]['positions'])
    comp = tsv(T, 't5_completeness.tsv')
    fr = S['floor_rationale']
    for r in recs(comp):
        c('T5', f'{r.gene} completeness', r.all_genomes, fr['gene_completeness'][r.gene], rel=1e-3)
        c('T5', f'{r.gene} worst position share', r.worst_position_call_rate, fr['worst_residue_position'][r.gene], rel=1e-3)
        for co in COH:
            c('T5', f'{r.gene} {co} completeness', getattr(r, co), fr['gene_completeness_by_cohort'][co][r.gene], rel=1e-3)
    mk = tf['markers']
    c('T5', 'allele columns exported', mk['allele_columns_exported'], S['exports']['alleles'])
    c('T5', 'residue columns exported', mk['residue_columns_exported'], S['exports']['residues'])
    c('T5', 'positions over 1%', mk['positions_losing_over_1pct'], S['residue_loss']['over_1pct'])
    c('T5', 'positions over 5%', mk['positions_losing_over_5pct'], S['residue_loss']['over_5pct'])
    # T6
    for r in recs(tsv(T, 't6_split_counts.tsv')):
        for k in ('tested', 'groups', 'columns', 'groups_compared', 'columns_compared'):
            c('T6', f'{r.cohort} {k}', getattr(r, k), S['pgroup_splits']['per_cohort'][r.cohort][k])
    rows = {x['p_group']: x for x in S['pgroup_mapping']['rows']}
    sg = tsv(T, 't6_split_groups.tsv')
    c('T6', 'P groups listed', len(sg), len(rows))
    for r in recs(sg):
        d = rows[r.p_group]
        c('T6', f'{r.p_group} tested', r.tested.replace('; ', ','), ','.join(x['allele'] for x in d['tested']))
        for a, f in zip(str(r.tested_freq_controls).split('; '), [x['freq'] for x in d['tested']]):
            c('T6', f'{r.p_group} freq', a, f, rel=1e-4)
        c('T6', f'{r.p_group} others', r.others_n, d['others_n'])
        c('T6', f'{r.p_group} others freq', r.others_freq, d['others_freq'])
        c('T6', f'{r.p_group} sum', r.group_freq_controls, d['ours_freq'])
        if d['compared']:
            c('T6', f'{r.p_group} ToMMo', r.tommo_freq, d['tommo_freq'])
    we = {x['allele']: x for x in S['worked_example']['rows']}
    for r in recs(tsv(T, 't6_dqa1.tsv')):
        d = we[r.p_group]
        for k, v in (('where', d['where']), ('copies_ours', d['copies_ours']), ('copies_tommo', d['copies_tommo']),
                     ('freq_ours', d['freq_ours']), ('freq_tommo', d['freq_tommo'])):
            c('T6', f'{r.p_group} {k}', getattr(r, k), v)
        c('T6', f'{r.p_group} members', '' if pd.isna(r.members) else r.members,
          ', '.join(f"{m['allele']} ({m['copies']})" for m in d['members']))


# ── association deck ────────────────────────────────────────────────────────
def association(c, A, S, RA, PG):
    s = tsv(A, 'a1_samples.tsv')
    for r in recs(s):
        d = S['samples'][r.comparison][r.cohort]
        for m in ('random', 'fixed'):
            for k in ('samples', 'cases', 'controls'):
                c('A1', f'{r.comparison} {r.cohort} {m} {k}', getattr(r, f'{m}_{k}'), d[m][k])
        c('A1', f'{r.comparison} {r.cohort} variance ratio', r.variance_ratio, S['null_models'][r.comparison][r.cohort])
        if r.comparison == 'all':
            for k in ('yes', 'no', 'unknown'):
                c('A1', f'{r.cohort} dvt {k}', getattr(r, f'dvt_{k}'), S['dvt'][r.cohort][k])
    g = tsv(A, 'a1_tests_by_gene.tsv')
    for co in COH:
        pg = {(x['gene'], x['class']): x['tests'] for x in S['per_gene'][co]}
        for r in recs(g[g.gene != 'All']):
            for cls in ('allele', 'residue'):
                c('A1', f'{co} {r.gene} {cls}', getattr(r, f'{cls}_{co}'), pg.get((r.gene, cls), 0))
    # A2
    th = {(x['case_set'], x['cohort'], x['class']): x for x in S['thresholds']}
    bt = {(x['case_set'], x['cohort'], x['class']): x for x in S['boot']}
    sc = tsv(A, 'a2_scans.tsv')
    c('A2', 'scans', len(sc), len(th))
    for r in recs(sc):
        k = (r.comparison, r.cohort, r['class'])
        for f in ('markers', 'mirror_pairs', 'tests', 'me_gec', 'p_gec', 'p_bonf', 'gec_blocks'):
            c('A2', f'{k} {f}', getattr(r, f), th[k][f])
        for f in ('draws', 'samples', 'cases', 'controls'):
            c('A2', f'{k} boot {f}', getattr(r, f'boot_{f}'), bt[k][f])
        for t in ('score', 'wald'):
            for f in ('p_maxt', 'p_maxt_{t}_lo', 'p_maxt_{t}_hi', 'fwer_gec', 'fwer_gec_{t}_lo', 'fwer_gec_{t}_hi',
                      'fwer_bonf', 'fwer_bonf_{t}_lo', 'fwer_bonf_{t}_hi', 'replicates_with_na'):
                key = f.format(t=t) if '{t}' in f else f'{f}_{t}'
                c('A2', f'{k} {key}', getattr(r, f'boot_{key}'), bt[k][key])
    ob = S['omnibus']
    for r in recs(tsv(A, 'a2_positions.tsv')):
        c('A2', f'{r.comparison} {r.cohort} positions tested', r.tested, ob['tested_by_scan'][f'{r.comparison}|residue|{r.cohort}'])
    for r in recs(tsv(A, 'a2_gec_check.tsv')):
        d = S['gec_check'][r['class']]
        for k, v in (('tests', d['tests']), ('blocks', d['blocks']), ('duplicates', d['duplicates']),
                     ('me_recomputed', d['me_reconstructed']), ('me_gec', d['me_gec']), ('max_block_diff', d['max_block_diff'])):
            c('A2', f'GEC check {r['class']} {k}', getattr(r, k), v)
    # A3
    a3 = tsv(A, 'a3_scans.tsv')
    for r in recs(a3[a3.comparison == 'all']):
        d = RA['scans'][r.test_type][r.cohort]
        c('A3', f'{r.cohort} {r.test_type} tests', r.tests, d['tests'])
        if r.test_type == 'position':
            c('A3', f'{r.cohort} position threshold', r.threshold_primary, d['p_bonf'])
            c('A3', f'{r.cohort} position sig', r.sig_primary, d['sig_bonf'])
        else:
            c('A3', f'{r.cohort} {r.test_type} gec', r.threshold_primary, d['p_gec'])
            c('A3', f'{r.cohort} {r.test_type} bonf', r.threshold_bonf, d['p_bonf'])
            c('A3', f'{r.cohort} {r.test_type} sig gec', r.sig_primary, d['sig_gec'])
            c('A3', f'{r.cohort} {r.test_type} sig bonf', r.sig_bonf, d['sig_bonf'])
    rand = {(x['case_set'], x['class'], x['cohort']): x for x in S['scans'] if x['model'] == 'random'}
    for r in recs(a3[a3.test_type != 'position']):
        d = rand[(r.comparison, r.test_type, r.cohort)]
        c('A3', f'{r.comparison} {r.cohort} {r.test_type} min p', r.min_p, d['min_p'])
        c('A3', f'{r.comparison} {r.cohort} {r.test_type} top', r.top, d['top_marker'])
        c('A3', f'{r.comparison} {r.cohort} {r.test_type} n gec', r.sig_primary, d['tests_gec'])
    for k, d in ob['min_p'].items():
        cs = k.split('|')[0]
        x = a3[(a3.comparison == cs) & (a3.test_type == 'position')]
        c('A3', f'{cs} smallest position P', x.min_p.min(), d['p'])
    top = tsv(A, 'a3_top.tsv')
    for k in ('allele', 'residue', 'position'):
        t = top[top.test_type == k].reset_index(drop=True)
        dk = RA['top'][k]
        c('A3', f'top {k} count', len(t), len(dk))
        for i, r in t.iterrows():
            d = dk[i]
            c('A3', f'top {k} {i + 1}', r.marker, d.get('marker', d.get('unit')))
            for co in COH:
                c('A3', f'top {k} {r.marker} p {co}', r[f'p_{co}'], d['p'][co])
                if k != 'position':
                    c('A3', f'top {k} {r.marker} or {co}', r[f'or_{co}'], d['or'][co])
            if k != 'position':
                c('A3', f'top {k} {r.marker} freq case', r.freq_case, d['freq_case'])
                c('A3', f'top {k} {r.marker} freq control', r.freq_control, d['freq_control'])
                c('A3', f'top {k} {r.marker} n3', r.n3_primary, d['gec'])
                c('A3', f'top {k} {r.marker} n3 bonf', r.n3_bonf, d['bonf'])
            else:
                c('A3', f'top position {r.marker} df', r.df, d['df'])
                c('A3', f'top position {r.marker} residues', r.residues, ', '.join(d['residues']))
                c('A3', f'top position {r.marker} n3', r.n3_primary, d['bonf'])
    # A4
    st = {(x['case_set'], x['class'], x['step']): x for x in S['steps'] if x['model'] == 'random'}
    for r in recs(tsv(A, 'a4_steps.tsv')):
        if r.test_type == 'position':
            continue
        d = st[(r.comparison, r.test_type, r.step)]
        c('A4', f'{r.comparison} {r.test_type} {r.step} lead', r.lead, d['lead'])
        c('A4', f'{r.comparison} {r.test_type} {r.step} entered', bool(r.entered), bool(d['entered']))
        for co in COH:
            c('A4', f'{r.comparison} {r.test_type} {r.step} p {co}', getattr(r, f'p_{co}'), d[f'p_{SHORT[co]}'])
            c('A4', f'{r.comparison} {r.test_type} {r.step} or {co}', getattr(r, f'or_{co}'), d[f'or_{SHORT[co]}'])
        c('A4', f'{r.comparison} {r.test_type} {r.step} n3', r.n3, d['gec'])
    af = tsv(A, 'a4_after.tsv')
    af = af[af.comparison == 'all']
    for co in COH:
        x = af[af.cohort == co].set_index('test_type')
        c('A4', f'{co} residue left', x.loc['residue', 'min_p'], RA['cond']['residue']['left'][co]['p'])
        c('A4', f'{co} residue left marker', x.loc['residue', 'top'], RA['cond']['residue']['left'][co]['marker'])
        c('A4', f'{co} position left', x.loc['position', 'min_p'], RA['cond']['position']['left'][co]['p'])
        c('A4', f'{co} position left unit', x.loc['position', 'top'], RA['cond']['position']['left'][co]['unit'])
        c('A4', f'{co} alleles given 57 V', x.loc['allele', 'min_p'], RA['summary']['alleles_given_57v'][co]['min_p'])
        c('A4', f'{co} alleles given 57 V sig', x.loc['allele', 'sig'], RA['summary']['alleles_given_57v'][co]['sig'])
    fa = json.loads((Path(A) / 'assoc_facts.json').read_text())
    for co in COH:
        for k in ('before', 'after'):
            c('A4', f'{co} B 114 {k}', fa['b114'][co][k], RA['cond']['position']['b114'][co][k])
    # A5
    gl = {x['lead_id']: x for x in S['gwas_leads']}
    for r in recs(tsv(A, 'a5_leads.tsv')):
        for co in COH:
            c('A5', f'{r.lead} p {co}', getattr(r, f'p_{co}'), gl[r.lead][f'p_{co}'])
            c('A5', f'{r.lead} or {co}', getattr(r, f'or_{co}'), gl[r.lead][f'or_{co}'])
            c('A5', f'{r.lead} freq {co}', getattr(r, f'freq_{co}'), gl[r.lead][f'a1_freq_{co}'])
    for r in recs(tsv(A, 'a5_link.tsv')):
        d = RA['gwas'][r.cohort]
        for k in ('lead_p', 'lead_given_57v', 'lead_given_pos57', 'tag', 'tag_r2', 'lead_given_tag', 'v57_p', 'v57_given_lead',
                  'hla_min_given_lead', 'hla_top_given_lead', 'most_explaining', 'r_57v'):
            c('A5', f'{r.cohort} {k}', getattr(r, k), d[k])
        d = RA['pos57_lead'][r.cohort]
        c('A5', f'{r.cohort} pos57', r.pos57_p, d['p'])
        c('A5', f'{r.cohort} pos57 given lead', r.pos57_given_lead, d['p_given_lead'])
        c('A5', f'{r.cohort} pos57 threshold', r.pos57_threshold, d['threshold'])
    lg = tsv(A, 'a5_lead_given.tsv')
    dl = [x for x in S['lead_given'] if x['model'] == 'random']
    c('A5', 'lead given rows', len(lg), len(dl))
    dk = {(x['case_set'], x['cohort'], x['conditioned_on']): x for x in dl}
    for r in recs(lg):
        d = dk[(r.comparison, r.cohort, r.conditioned_on)]
        for k in ('lead_p_marginal', 'lead_p_given', 'lead_or_given', 'marker_p_marginal', 'marker_p_given_lead'):
            c('A5', f'{r.comparison} {r.cohort} {r.conditioned_on} {k}', getattr(r, k), d[k])
    # A6
    lam = tsv(A, 'a6_lambda.tsv')
    key = {'GWAS lead': 'lead', 'DQB1 57 V': 'v57', 'position 57': 'pos57'}
    for r in recs(lam):
        if r.test_type in ('allele', 'residue', 'position') and r.conditioning in key:
            d = RA['lambda'][key[r.conditioning]][r.test_type][r.cohort]
            c('A6', f'{r.test_type} {r.conditioning} {r.cohort} tests', r.tests, d['tests'])
            c('A6', f'{r.test_type} {r.conditioning} {r.cohort} before', r.before, d['before'])
            c('A6', f'{r.test_type} {r.conditioning} {r.cohort} after', r.after, d['after'])
        elif r.test_type in ('allele', 'residue', 'position') and r.conditioning == 'none':
            d = RA['lambda_marginal'][r.test_type][r.cohort]
            c('A6', f'{r.test_type} {r.cohort} alone', r.before, d['lambda'])
            c('A6', f'{r.test_type} {r.cohort} tests', r.tests, d['tests'])
        elif r.test_type.startswith('position,'):
            d = RA['lambda_position_df'][r.cohort]['1' if '1 df' in r.test_type else '2+']
            c('A6', f'{r.test_type} {r.cohort} n', r.tests, d['n'])
            c('A6', f'{r.test_type} {r.cohort} before', r.before, d['before'])
            c('A6', f'{r.test_type} {r.cohort} given', r.after, d['given_v57'])
        elif r.test_type == 'GWAS SNPs in the MHC':
            d = RA['lambda_marginal']['gwas'][r.cohort]
            c('A6', f'GWAS MHC {r.cohort}', r.before, d['lambda'], rel=1e-5)
            c('A6', f'GWAS MHC {r.cohort} SNPs', r.tests, d['tests'])
    # A7
    for r in recs(tsv(A, 'a7_v57.tsv')):
        d = RA['summary'][r.cohort]
        for k in ('or', 'l95', 'u95', 'p', 'or_fixed', 'l95_fixed', 'u95_fixed', 'p_fixed', 'gec'):
            c('A7', f'{r.cohort} {k}', getattr(r, k), d[k])
        o = RA['pos57']['omnibus'][r.cohort]
        c('A7', f'{r.cohort} pos57 chisq', r.chisq_pos57, o['chisq'])
        c('A7', f'{r.cohort} pos57 p', r.p_pos57, o['p'])
        c('A7', f'{r.cohort} pos57 threshold', r.threshold_pos57, o['threshold'])
    for r in recs(tsv(A, 'a7_genotypes.tsv')):
        d = RA['pos57']['genotypes'][r.cohort]['case' if r.group == 'cases' else 'control']
        for k, v in (('analyzed', d['n']), ('called', d['called']), ('undetermined', d['missing']), ('n0', d['n0']),
                     ('n1', d['n1']), ('n2', d['n2'])):
            c('A7', f'{r.cohort} {r.group} {k}', getattr(r, k), v)
    rs = tsv(A, 'a7_residues.tsv')
    for r in recs(rs[rs.level == 'residue']):
        d = RA['pos57']['residues'][r.cohort][r.residue]
        for k, v in (('n_alleles', d['n_alleles']), ('copies_case', d['case']['copies']), ('copies_control', d['control']['copies']),
                     ('or', d['or']), ('p', d['p'])):
            c('A7', f'{r.cohort} {r.residue} {k}', getattr(r, k), v)
    for r in recs(rs[(rs.level == 'allele') & rs.name.str.contains(r'\*', regex=True)]):
        d = next(x for x in RA['pos57']['alleles'][r.cohort][r.residue] if x['name'] == r.name)
        c('A7', f'{r.cohort} {r.name} copies case', r.copies_case, d['case']['copies'])
        c('A7', f'{r.cohort} {r.name} copies control', r.copies_control, d['control']['copies'])
        if d.get('tested'):
            c('A7', f'{r.cohort} {r.name} p', r.p, d['p'])
    # A8
    for r in recs(tsv(A, 'a8_strata.tsv')):
        d = RA['supp']['marker']['DQB1:57:V'][r.comparison][r.cohort]
        for k in ('freq_case', 'freq_control', 'or', 'l95', 'u95', 'p'):
            c('A8', f'{r.comparison} {r.cohort} {k}', getattr(r, k), d[k])
        n = RA['supp']['strata'][r.comparison][r.cohort]
        c('A8', f'{r.comparison} {r.cohort} n case', r.n_case, n['n_case'])
        c('A8', f'{r.comparison} {r.cohort} n control', r.n_control, n['n_control'])
    for r in recs(a3[a3.comparison != 'all']):
        d = RA['supp']['scan'][r.comparison][r.cohort]
        c('A8', f'{r.comparison} {r.cohort} {r.test_type} sig', r.sig_primary, d[r.test_type])
    for r in recs(tsv(A, 'a8_nodvt_positions.tsv')):
        for co in COH:
            c('A8', f'{r.unit} {co}', getattr(r, f'p_{co}'), RA['supp']['pos_all3']['p'][r.unit][co])
    for r in recs(tsv(A, 'a8_lrt_vs_score.tsv')):
        d = RA['supp']['lrt_vs_score']['cohorts'][r.cohort]
        c('A8', f'{r.cohort} LRT', r.p_lrt_fixed, d['p_lrt'])
        c('A8', f'{r.cohort} score', r.p_score_random, d['p_score'])
    nd = {x['kind']: x for x in RA['supp']['nodvt_dqb1']['rows']}
    for r in recs(tsv(A, 'a8_nodvt_dqb1.tsv')):
        d = nd[r.kind]
        c('A8', f'{r.kind} n', r.n_alleles, d['n_alleles'])
        c('A8', f'{r.kind} freq case', r.freq_case, d['freq_case'])
        c('A8', f'{r.kind} freq control', r.freq_control, d['freq_control'])
        for q in ('57', '87', '125'):
            c('A8', f'{r.kind} residues {q}', '' if pd.isna(getattr(r, f'res_{q}')) else getattr(r, f'res_{q}'),
              ', '.join(d['residues'][q]))
    ld = fa['ld_57_87_without_dvt']
    for k in ('r2', 'r2_ceiling', 'dprime'):
        c('A8', f'LD {k}', ld[k], RA['supp']['nodvt_dqb1']['ld'][k])
    g57 = tsv(A, 'a8_given57v.tsv').set_index('cohort')
    for co in COH:
        d = RA['supp']['given57v'][co]
        for k in ('allele', 'residue', 'position'):
            c('A8', f'{co} given57v {k} before', g57.loc[co, f'{k}_before'], d[k]['before'])
            c('A8', f'{co} given57v {k} after', g57.loc[co, f'{k}_after'], d[k]['after'])
            c('A8', f'{co} given57v {k} top', g57.loc[co, f'{k}_top'], d[k]['top'])
            c('A8', f'{co} given57v {k} top p', g57.loc[co, f'{k}_min_p'], d[k]['top_p'])
        for m in ('DQB1:87:Y', 'DPB1*05:01', 'DQB1*05:03', 'DQB1*05:03|all', 'DQB1:125'):
            for k in ('p', 'p_given') + (('or', 'or_given') if m != 'DQB1:125' else ()):
                c('A8', f'{co} {m} {k}', g57.loc[co, f'{m}|{k}'], d[m][k])
        for m, k in (('DQB1:87:Y', 'p_given_87y'), ('DQB1:87:Y', 'or_given_87y')):
            pass
    for r in recs(tsv(A, 'a8_57v_given_87y.tsv')):
        d = RA['supp']['cond'][r.cohort]['DQB1:87:Y']
        c('A8', f'{r.cohort} 57V marginal', r.p, d['p_marginal'])
        c('A8', f'{r.cohort} 57V given 87Y', r.p_given_87y, d['p_given'])
        c('A8', f'{r.cohort} 57V OR given 87Y', r.or_given_87y, d['or_given'])
    for k, d in fa['spa_calls'].items():
        c('A8', f'{k} SPA calls', d['calls'], S['spa'][k]['calls_na'])
        c('A8', f'{k} SPA calls spa', d['calls_spa'], S['spa'][k]['calls_spa'])
    # A9
    pub = {(x['study'], x['allele']): x for x in S['replication']['published']}
    for r in recs(tsv(A, 'a9_published.tsv')):
        d = pub[(r.study.split()[0].lower(), r.allele)]
        for k in ('case_car', 'case_n', 'ctrl_car', 'ctrl_n', 'or_', 'l95', 'u95', 'p', 'pc'):
            c('A9', f'{r.study} {r.allele} {k}', getattr(r, k), d[k])
    ours = {(x['cohort'], x['set'], x['allele']): x for x in S['replication']['ours']}
    for r in recs(tsv(A, 'a9_carriers.tsv')):
        d = ours[(r.cohort, r.comparison, r.allele)]
        for k, v in (('case_car', d['case_car']), ('case_n', d['case_n']), ('ctrl_car', d['ctrl_car']), ('ctrl_n', d['ctrl_n']),
                     ('or', d['or_']), ('l95', d['l95']), ('u95', d['u95']), ('p', d['p']),
                     ('power_published', d['power_published']), ('detectable_or80', d['detectable_or80'])):
            c('A9', f'{r.cohort} {r.comparison} {r.allele} {k}', getattr(r, k), v)
    for a, d in fa['pop_reported'].items():
        c('A9', f'{a} ToMMo', d['tommo'], S['replication']['population'][a]['allele_freq_tommo'])
        c('A9', f'{a} controls', d['controls_pgroup'], S['replication']['population'][a]['allele_freq_ours_controls'])
    # A10
    cc = {(x['case_set'], x['cohort'], x['class']): x for x in S['concordance']}
    for r in recs(tsv(A, 'a10_concordance.tsv')):
        d = cc[(r.comparison, r.cohort, r.test_type)]
        for k in ('tests', 'r_nlp', 'r_beta', 'median_ratio'):
            c('A10', f'{r.comparison} {r.cohort} {r.test_type} {k}', getattr(r, k), d[k])
    om = tsv(A, 'a10_one_model_only.tsv')
    c('A10', 'one model only rows', len(om), len(S['one_model_only']))
    for r, d in zip(recs(om), S['one_model_only']):
        c('A10', f'one model {r.marker}', r.marker, d['marker'])
        c('A10', f'one model {r.marker} in', r.significant_in, d['sig_in'])
        c('A10', f'one model {r.marker} p fixed', r.p_fixed, d['p_fixed'])
        c('A10', f'one model {r.marker} p random', r.p_random, d['p_random'])
    for k, v in S['mirror_summary'].items():
        c('A10', f'mirror {k}', fa['mirror_full'][k], v)
    pl = tsv(A, 'a10_platform.tsv')
    for r in recs(pl):
        bg = S['platform_by_group'][r.cohort]
        c('A10', f'{r.cohort} {r.platform} controls', r.controls, bg['0'].get(r.platform, 0))
        c('A10', f'{r.cohort} {r.platform} cases', r.cases_all, bg['1'].get(r.platform, 0))
        t = S['platform_within_cases'][r.cohort]['table'].get(r.platform)
        if t:
            for k in ('yes', 'no', 'unknown'):
                c('A10', f'{r.cohort} {r.platform} dvt {k}', getattr(r, f'cases_dvt_{k}'), t[k])
    for r in recs(tsv(A, 'a10_platform_chi2.tsv')):
        d = S['platform_within_cases'][r.cohort]
        for k in ('chi2', 'df', 'p'):
            c('A10', f'{r.cohort} platform {k}', getattr(r, k), d[k])
    pp = tsv(A, 'a10_pos57_pgroups.tsv').set_index('p_group')
    c('A10', 'P groups at 57', len(pp), len(PG))
    for r in recs(PG):
        x = pp.loc[r["name"]]
        c('A10', f'{r.p_group} residue', x.residue_57, r.residue_57)
        c('A10', f'{r.p_group} cases', x.copies_case, r.copies_case)
        c('A10', f'{r.p_group} controls', x.copies_control, r.copies_control)
        c('A10', f'{r.p_group} ToMMo', x.tommo, r.tommo)
        c('A10', f'{r.p_group} 1KG', x.kg, r.kg_jpt)


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for k in ('typing', 'assoc', 'typing_stats', 'assoc_stats', 'assoc_results', 'assoc_pgroups', 'out'):
        ap.add_argument('--' + k.replace('_', '-'), required=True)
    a = ap.parse_args()
    c = Cmp()
    typing(c, a.typing, json.loads(Path(a.typing_stats).read_text()))
    association(c, a.assoc, json.loads(Path(a.assoc_stats).read_text()), json.loads(Path(a.assoc_results).read_text()),
                pd.read_csv(a.assoc_pgroups, sep='\t', dtype=str).fillna('').pipe(
                    lambda d: d.astype({k: int for k in ('copies_case', 'copies_control', 'tommo', 'kg_jpt')})))
    sheets = {'T1': 'T1 Typing run', 'T2': 'T2 DRB3, DRB4, DRB5', 'T3': 'T3 Completeness', 'T4': 'T4 ToMMo comparison',
              'T5': 'T5 Marker QC', 'T6': 'T6 P groups', 'A1': 'A1 Design', 'A2': 'A2 Thresholds',
              'A3': 'A3 Q1 Association', 'A4': 'A4 Q2 Conditioning', 'A5': 'A5 Q3 GWAS lead', 'A6': 'A6 Inflation',
              'A7': 'A7 Summary', 'A8': 'A8 DVT', 'A9': 'A9 Reported alleles', 'A10': 'A10 Appendix'}
    checks = []
    for k, name in sheets.items():
        if not c.n[k]:
            raise SystemExit(f'deck_check: nothing compared for {name}')
        deck = 'HLA typing deck' if k.startswith('T') else 'HLA association deck'
        checks.append({'id': f'deck|{k}', 'group': 'Deck agreement', 'sheet': name,
                       'what': f'every quantity the sheet shares with the {deck}, against the data the deck was built from',
                       'route': 'deck data files (below)', 'compared': c.n[k], 'mismatches': len(c.bad[k]),
                       'passed': not c.bad[k], 'detail': c.bad[k][:50]})
    files = [{'deck': 'HLA typing', 'path': str(Path(a.typing_stats).resolve()), 'md5': md5(a.typing_stats)}] + \
            [{'deck': 'HLA association', 'path': str(Path(p).resolve()), 'md5': md5(p)}
             for p in (a.assoc_stats, a.assoc_results, a.assoc_pgroups)]
    out = {'all_passed': all(x['passed'] for x in checks), 'checks': checks, 'deck_files': files}
    Path(a.out).write_text(json.dumps(out, indent=1, ensure_ascii=False) + '\n')
    n, m = sum(x['compared'] for x in checks), sum(x['mismatches'] for x in checks)
    print(f'[deck_check] {n} values compared, {m} mismatches')
    for x in checks:
        for d in x['detail'][:8]:
            print(f"  {x['sheet']}: {d['item']}: workbook {d['workbook']} | deck {d['deck']}")


if __name__ == '__main__':
    main()
