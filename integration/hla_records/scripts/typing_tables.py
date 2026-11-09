#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The Typing part of the workbook (sheets T1-T8), assembled from the
#           finished hla.typing.v2 run. The definitions are the typing deck's
#           (its collect_stats.py), so every number shared with a slide is the
#           same quantity; the rules (gene sets, floors, P groups, releases) come
#           from the typing component's own config.py / nomenclature.py.
#             T1 typing run          run facts; samples by group x platform
#             T2 DRB3/4/5            copy number against DRB1
#             T3 completeness        per gene: states, basis, reads, confidence, tier
#             T4 ToMMo comparison    harmonization; per gene and panel; DRB1*15:13;
#                                    homozygosity, confidence vs. panel
#             T5 marker QC           encoding; stepwise record; per gene; floors
#             T6 P groups            tested columns sharing a P group; DQA1 example
#             T7 alleles (detail)    every exported allele column
#             T8 residues (detail)   every exported residue column
#           Nothing per sample is written: every table is a count or a share.
# Inputs  : --typing (hla.typing.v2 dir), --model-inputs (wgs.auto.par 12_model_inputs)
# Outputs : --out/t*.tsv, --out/typing_facts.json
# Component: hla_records
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hr_common import COHORTS, die, log, write_json, write_tsv  # noqa: E402

STAGES = [('typed', 'exported', 'Typed by HLA-HD', 'columns that vary in the cohort'),
          ('genes', 'genes', 'DRB3/4/5 excluded', 'their presence follows DRB1'),
          ('callrate', 'call_rate', 'Call rate ≥ 95% in cases and in controls', 'share of samples with a genotype'),
          ('maf', 'maf', 'MAF ≥ 0.5%', 'copies of the less common form')]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--typing', required=True, help='hla.typing.v2 directory')
    p.add_argument('--model-inputs', required=True, help='wgs.auto.par/results/12_model_inputs')
    p.add_argument('--out', required=True)
    return p.parse_args()


class Typing:
    """Every typing table read once."""

    def __init__(self, a):
        self.root = Path(a.typing)
        sys.path.insert(0, str(self.root / 'scripts'))
        import config as C            # noqa: E402  (the typing component's rules)
        import nomenclature as N      # noqa: E402
        self.C, self.N = C, N
        R = self.R = self.root / 'results'
        self.mi = Path(a.model_inputs)
        self.prov = json.loads((R / '00.inputs/provenance.json').read_text())
        self.si = pd.read_csv(R / '00.inputs/sample_inputs.tsv', sep='\t', dtype={'sample_id': str})
        self.sq = pd.read_csv(R / '02.evidence/sample_qc.tsv', sep='\t', dtype={'sample_id': str})
        self.ev = pd.read_csv(R / '02.evidence/call_evidence.tsv', sep='\t', dtype=str).fillna('')
        self.gs = pd.read_csv(R / '02.evidence/gene_summary.tsv', sep='\t')
        self.drb = pd.read_csv(R / '02.evidence/drb_consistency.tsv', sep='\t', dtype=str)
        self.rec = pd.read_csv(R / '01.calls/drb_reconciliation.tsv', sep='\t', dtype=str).fillna('')
        self.calls = pd.read_csv(R / '01.calls/calls_long.tsv', sep='\t', dtype=str).fillna('')
        self.loc = pd.read_csv(R / '03.panels/locus_concordance.tsv', sep='\t')
        self.ac = pd.read_csv(R / '03.panels/allele_concordance.tsv', sep='\t')
        self.hz = pd.read_csv(R / '03.panels/homozygosity.tsv', sep='\t')
        self.cv = pd.read_csv(R / '03.panels/confidence_vs_panel.tsv', sep='\t')
        self.unc = pd.read_csv(R / '03.panels/unconfirmed_alleles.tsv', sep='\t')
        self.shift = pd.read_csv(R / '03.panels/shift_direction.tsv', sep='\t')
        self.tiers = pd.read_csv(R / '04.confidence/locus_tiers.tsv', sep='\t')
        self.ad = pd.read_csv(R / '05.exports/allele_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
        self.rd = pd.read_csv(R / '05.exports/residue_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
        self.pgmap = pd.read_csv(R / '05.exports/allele_pgroup_map.tsv', sep='\t')
        self.match = pd.read_csv(R / '05.exports/allele_match_map.tsv', sep='\t')
        self.sites = pd.read_csv(R / '05.exports/residue_sites.tsv', sep='\t')
        self.A = R / '06.assoc_inputs'
        self.mk = pd.read_csv(self.A / 'markers.tsv', sep='\t', dtype={'position': str, 'residue': str})
        self.qc = pd.read_csv(self.A / 'qc_record.tsv', sep='\t')
        self.summ = json.loads((self.A / 'summary.json').read_text())
        self.pheno = {c: pd.read_csv(self.mi / c / 'phenotype' / 'pheno.tsv', sep='\t', dtype={'IID': str})
                      for c in COHORTS}

    def controls_full(self):
        ph = self.pheno['full_mainland']
        return ph.loc[ph.PHENO1 == 0, 'IID']


# ── T1 ───────────────────────────────────────────────────────────────────────
def t1(T):
    C, si, sq, ev = T.C, T.si, T.sq, T.ev
    run = [('Genomes typed', len(si), 'cohort.n'),
           ('Controls (AGP3K)', int((si.group == C.CONTROL_GROUP).sum()), 'cohort.groups.AGP3K'),
           ('Cases (CTEPH)', int((si.group != C.CONTROL_GROUP).sum()), 'cohort.groups.PH'),
           ('Runs validated', int((si.validation == 'ok').sum()), 'cohort.validation_ok'),
           ('HLA-HD version', T.prov['hlahd_version'], 'run.hlahd_version'),
           ('IPD-IMGT/HLA release of the dictionary', '3.64.0', 'run.imgt'),
           ('Gene split file release', '3.50.0', 'run.gene_split'),
           ('Read regions extracted', T.prov['regions'], 'run.regions'),
           ('Genes typed by HLA-HD', int(T.calls.gene.nunique()), 'run.genes_typed'),
           ('Genes reported', len(C.REPORTED_GENES), 'run.genes_reported'),
           ('Genes tested in the association', len(C.ASSOC_GENES), 'assoc_basis.genes'),
           ('HLA read pairs per sample, median', int(si.read_pairs.median()), 'cohort.read_pairs_median'),
           ('HLA read pairs, 1st percentile (outlier limit)', int(np.percentile(si.read_pairs, C.OUTLIER_READ_PAIRS_PCTL)),
            'cohort.read_pairs_p1'),
           ('HLA read pairs, minimum', int(si.read_pairs.min()), 'cohort.read_pairs_min'),
           ('HLA read pairs, maximum', int(si.read_pairs.max()), 'cohort.read_pairs_max'),
           ('Calls (genomes × reported genes)', len(ev), 'calls_reported'),
           ('Sample outliers (listed, never removed)', int(sq.outlier.sum()), 'samples.outliers')]
    run = pd.DataFrame(run, columns=['item', 'value', 'deck_key'])
    # per group x platform: the platform groups are never merged (WGS_Platform as recorded)
    two = ev.state.isin(['heterozygous', 'homozygous', 'present', 'absent'])
    e = ev.assign(two=two, high=ev.confidence == 'high')
    per = e.groupby('sample_id').agg(calls=('gene', 'size'), two=('two', 'sum'), high=('high', 'sum'))
    q = sq.set_index('sample_id').join(per)
    q['reasons'] = q.outlier_reason.fillna('')
    rows = []

    def summarise(d, group, platform):
        rr = d.reasons.str.split(';')
        return {'group': group, 'platform': platform, 'samples': len(d),
                'male': int((d.sex == 'M').sum()), 'female': int((d.sex == 'F').sum()),
                'read_pairs_median': float(d.read_pairs.median()), 'read_pairs_min': int(d.read_pairs.min()),
                'depth_median': float(d.observed_depth.median()), 'depth_min': float(d.observed_depth.min()),
                'depth_max': float(d.observed_depth.max()),
                'outliers': int(d.outlier.sum()),
                'out_read_pairs': int(rr.map(lambda x: 'read_pairs<p1' in x).sum()),
                'out_no2field': int(rr.map(lambda x: 'no_2field>=3' in x).sum()),
                'out_failed_obligate': int(rr.map(lambda x: 'failed_obligate' in x).sum()),
                'calls': int(d.calls.sum()), 'genotypes': int(d.two.sum()), 'genotype_share': d.two.sum() / d.calls.sum(),
                'high_conf': int(d.high.sum()), 'high_share': d.high.sum() / d.calls.sum()}
    ctrl = q[q.group == C.CONTROL_GROUP]
    for pf, d in ctrl.groupby('platform'):
        rows.append(summarise(d, 'Controls', pf))
    case = q[q.group != C.CONTROL_GROUP]
    order = ['DNBSeq-T7 30x', 'NovaSeq 30x', 'DNBSeq-G400RS 30x', 'DNBseq-G400RS 15x']
    for pf in order + sorted(set(case.platform) - set(order)):
        d = case[case.platform == pf]
        if len(d):
            rows.append(summarise(d, 'Cases', pf))
    rows.append(summarise(case, 'Cases', 'All cases'))
    rows.append(summarise(q, 'All', 'All genomes'))
    plat = pd.DataFrame(rows)
    if plat[plat.platform == 'All genomes'].samples.iloc[0] != len(si) or len(per) != len(si):
        die('T1: samples do not add up')
    if int(plat[plat.platform == 'All genomes'].outliers.iloc[0]) != int(sq.outlier.sum()):
        die('T1: outliers do not add up')
    reasons = sq.outlier_reason[sq.outlier].str.split(';').explode().value_counts().to_dict()
    return run, plat, reasons


# ── T2 ───────────────────────────────────────────────────────────────────────
def t2(T):
    C, drb, rec, ev, calls, ad = T.C, T.drb.copy(), T.rec, T.ev, T.calls, T.ad
    drb['samples'] = drb.samples.astype(int)
    rows = []
    for g in C.PARALOGUES:
        d = drb[(drb.gene == g) & (drb.expected != 'unknown')]
        flags = rec[rec.gene == g].flag.value_counts().to_dict()
        tot = ad[[c for c in ad.columns if c.split('*')[0] == g]].sum(axis=1, min_count=1)
        r = rec[rec.gene == g].set_index('sample_id')
        e = ev[ev.gene == g].set_index('sample_id')
        exp = pd.to_numeric(r.expected, errors='coerce')
        ids = e.index[e.state.isin(['present', 'absent']) & exp.reindex(e.index).notna()]
        rows.append({'gene': g, 'with_drb1': int(d.samples.sum()),
                     'agree': int(d[d.expected == d.observed].samples.sum()),
                     'two_where_one': int(d[(d.expected == '1') & (d.observed == '2')].samples.sum()),
                     'present_where_none': int(d[(d.expected == '0') & (d.observed != '0')].samples.sum()),
                     'no_flag': int(flags.get('', 0)), 'reconciled': int(flags.get('drb_reconciled', 0)),
                     'unexpected': int(flags.get('drb_unexpected', 0)), 'missing': int(flags.get('drb_missing', 0)),
                     'drb1_unknown': int(flags.get('drb1_unknown', 0)),
                     'with_result': len(ids), 'exported_match': int((tot[ids] == exp[ids]).sum()),
                     'exported_blank': int(tot.isna().sum())})
    out = pd.DataFrame(rows)
    out['agree_share'] = out.agree / out.with_drb1
    out['no_flag_share'] = out.no_flag / out.with_drb1
    if (out.exported_match != out.with_result).any():
        die('T2: exported DRB copies disagree with DRB1')
    # raw HLA-HD calls of two alleles where DRB1 predicts one or two copies, with >= 10 unique reads
    p = calls[calls.gene.isin(C.PARALOGUES) & (calls.raw == 'two')].copy()
    u1, u2 = pd.to_numeric(p.unique_1), pd.to_numeric(p.unique_2)
    keep = (u1 + u2) >= C.MIN_UNIQUE_READS
    p, ms = p[keep], (np.minimum(u1, u2) / (u1 + u2))[keep]
    two = pd.DataFrame([{'expected_copies': int(e), 'two_allele_calls': int((p.drb_expected == e).sum()),
                         'minor_below_15pct': int(((p.drb_expected == e) & (ms < C.MIN_MINOR_SHARE)).sum())}
                        for e in ('1', '2')])
    return out, two


# ── T3 ───────────────────────────────────────────────────────────────────────
def imgt_alleles(T):
    import csv
    from collections import Counter
    cnt = Counter()
    with open(T.C.IMGT_HISTORY) as fh:
        rd = csv.reader(ln for ln in fh if not ln.startswith('#'))
        next(rd)
        for r in rd:
            nm = r[1].strip()
            if '*' in nm and nm != 'NA':
                cnt[nm.split('*')[0]] += 1
    return cnt


def t3(T):
    C, ev, gs, tiers = T.C, T.ev, T.gs.set_index('gene'), T.tiers.set_index('gene')
    cnt = imgt_alleles(T)
    two = ev.state.isin(['heterozygous', 'homozygous', 'present', 'absent'])
    res = ~ev.state.isin(C.NO_CALL_STATES)
    rows = []
    for g in C.REPORTED_GENES:
        r, t = gs.loc[g], tiers.loc[g]
        m = ev.gene == g
        rows.append({'gene': g, 'class': r['class'], 'assoc': g in C.ASSOC_GENES, 'imgt_alleles': int(cnt[g]),
                     'calls': int(m.sum()), 'heterozygous': int(r.state_heterozygous),
                     'homozygous': int(r.state_homozygous), 'present': int(r.state_present),
                     'absent': int(r.state_absent), 'ambiguous': int(r.state_ambiguous),
                     'no_result': int((m & ~res).sum()), 'unresolved': int(r.state_unresolved),
                     'failed': int(r.state_failed), 'missing': int(r.state_missing),
                     'genotypes': int((m & two).sum()), 'genotype_share': float(two[m].mean()),
                     'result_share': float(res[m].mean()),
                     'basis_agreed': int(r.basis_agreed), 'basis_prior': int(r.basis_prior),
                     'basis_pgroup': int(r.basis_pgroup), 'weak_prior': int(r.weak_prior),
                     'median_min_reads': float(r.median_min_reads), 'low_support': int(r.low_support),
                     'het_assessed': int(r.het_assessed), 'imbalanced': int(r.imbalanced),
                     'conf_high': int(r.conf_high), 'conf_medium': int(r.conf_medium), 'conf_low': int(r.conf_low),
                     'conf_no_call': int(r.conf_no_call), 'share_high': float(r.share_high),
                     'tier': int(t.tier), 'criteria_failed': '' if pd.isna(t.criteria_failed) else str(t.criteria_failed),
                     'external_panel': '' if pd.isna(t.external_panel) else str(t.external_panel)})
    out = pd.DataFrame(rows)
    # the deck's totals over allele slots: the second slot of a one-copy call (DRB3/4/5 present) is not an allele
    al = pd.concat([ev[['gene', 'basis_1', 'allele_1', 'prior_ratio_1', 'copies']].set_axis(
                        ['gene', 'basis', 'allele', 'q', 'copies'], axis=1),
                    ev[['gene', 'basis_2', 'allele_2', 'prior_ratio_2', 'copies']].set_axis(
                        ['gene', 'basis', 'allele', 'q', 'copies'], axis=1).assign(second=1)])
    al = al[~((al.second == 1) & (al.copies == '1'))]
    al = al[al.basis != '']
    q = pd.to_numeric(al.q, errors='coerce')
    h = ev[(ev.balance_assessed == 'True') & ev.gene.isin(C.OBLIGATE_GENES)]
    totals = {'calls': len(ev), 'result': int(res.sum()), 'genotypes': int(two.sum()),
              'slots': len(al), 'agreed': int((al.basis == 'agreed').sum()), 'prior': int((al.basis == 'prior').sum()),
              'prior_decisive': int(((al.basis == 'prior') & (q >= C.PRIOR_DOMINANCE)).sum()),
              'pgroup': int((al.basis == 'pgroup').sum()), 'none': int((al.basis == 'none').sum()),
              'het_assessed_obligate': len(h), 'imbalanced_obligate': int((h.imbalanced == 'True').sum()),
              'imgt_alleles': int(sum(cnt.values())),
              'states': ev.state.value_counts().to_dict(), 'confidence': ev.confidence.value_counts().to_dict(),
              'tier_counts': T.tiers.tier.value_counts().sort_index().to_dict()}
    if out.calls.sum() != totals['calls'] or out.genotypes.sum() != totals['genotypes']:
        die('T3: per-gene calls do not add up')
    return out, totals


# ── T4 ───────────────────────────────────────────────────────────────────────
def coverage_and_mapping(T):
    """Presence of the controls' alleles in ToMMo and the P-group mapping per compared gene (the typing
    deck's panel_coverage and pgroup_mapping)."""
    C, N = T.C, T.N
    loc = T.loc[T.loc.panel == 'ToMMo'].set_index('gene')
    vint = int(loc.panel_vintage.iloc[0])
    a = T.ac[(T.ac.panel == 'ToMMo') & (T.ac.count_ours > 0)].copy()
    rel = a.allele.map(N.release_of)
    a['state'] = np.where(a.count_panel > 0, 'in_panel', np.where(rel.isna() | (rel > vint), 'too_new', 'absent'))
    X = T.ad.loc[T.controls_full()]
    seen = X.sum()
    tm = pd.read_csv(C.TOMMO)
    acx = T.ac[T.ac.panel == 'ToMMo'].set_index(['gene', 'allele'])
    rows = []
    for g in [g for g in loc.index if int(loc.loc[g, 'frequency_comparable']) == 1]:
        d = a[a.gene == g]
        alle = d.groupby('state').allele.count()
        cop = d.groupby('state').count_ours.sum()
        cols = [c for c in X.columns if c.split('*')[0] == g and not c.endswith('N')]
        s = X[cols].sum()
        s = s[s > 0]
        if int(s.sum()) != int(loc.loc[g, 'n_chr_ours']):
            die(f'T4 {g}: controls copies differ from the comparison denominator')
        keys = pd.Series(s.index.map(N.pg_of), index=s.index)
        size = keys.map(keys.value_counts())
        t_ = tm[tm.gene == g]
        nonnull = ~t_.allele.str.rstrip('P').str.endswith('N')
        cmp_ = acx.loc[g]
        both = cmp_[(cmp_.count_ours > 0) & (cmp_.count_panel > 0)]
        row = {'gene': g, 'two_field': int(((seen > 0) & (seen.index.str.split('*').str[0] == g)
                                            & ~seen.index.str.endswith('N')).sum()),
               'ours_keys': int(keys.nunique()), 'ours_shared': int((size > 1).sum()),
               'shared_groups': int(keys[size > 1].nunique()),
               'tommo_p': int((t_.allele.str.endswith('P') & nonnull).sum()),
               'tommo_plain': int((~t_.allele.str.endswith('P') & nonnull).sum()),
               'tommo_keys': int((cmp_.count_panel > 0).sum()), 'both': int(len(both)), 'either': int(len(cmp_)),
               'r_shared': float(np.corrcoef(both.freq_ours, both.freq_panel)[0, 1]) if len(both) > 2 else np.nan}
        for k in ('in_panel', 'absent', 'too_new'):
            row[f'alleles_{k}'] = int(alle.get(k, 0))
            row[f'copies_{k}'] = int(cop.get(k, 0))
        if row['ours_keys'] != row['alleles_in_panel'] + row['alleles_absent'] + row['alleles_too_new']:
            die(f'T4 {g}: P-group keys do not add up')
        if row['either'] != int(loc.loc[g, 'n_alleles']) or row['both'] != int(loc.loc[g, 'n_shared']):
            die(f'T4 {g}: either/both differ from the comparison')
        rows.append(row)
    return pd.DataFrame(rows), vint


def t4(T):
    C, N = T.C, T.N
    cov, vint = coverage_and_mapping(T)
    loc = T.loc.copy()
    keep = ['panel', 'gene', 'frequency_comparable', 'note', 'n_samples', 'n_chr_ours', 'n_chr_panel', 'n_alleles',
            'n_shared', 'r', 'D', 'excess_D', 'unconfirmed', 'unconfirmed_assessable', 'unconfirmed_after_vintage',
            'null_chr_ours', 'null_chr_panel']
    genes = loc[keep].merge(cov, on='gene', how='left')
    genes.loc[genes.panel != 'ToMMo', [c for c in cov.columns if c != 'gene']] = np.nan
    genes['panel'] = genes.panel.replace({'ToMMo': 'ToMMo 61KJPN'})
    # harmonization (slide 15)
    gl = list(cov.gene)
    lt = T.loc[(T.loc.panel == 'ToMMo') & T.loc.gene.isin(gl)]
    tm = pd.read_csv(C.TOMMO)
    tm = tm[tm.gene.isin(gl)]
    null = tm.allele.str.rstrip('P').str.endswith('N')
    named_p = tm.allele.str.endswith('P') & ~null
    seen = T.ad.loc[T.controls_full()].sum()
    seen = seen[(seen > 0) & seen.index.str.split('*').str[0].isin(gl) & ~seen.index.str.endswith('N')]
    merged = {a: N.pg_of(a) for a in seen.index if N.pg_of(a) != a}
    ex_merge = max(merged, key=lambda a: seen[a])
    acb = T.ac[(T.ac.panel == 'ToMMo') & T.ac.gene.isin(gl) & (T.ac.count_ours > 0) & (T.ac.count_panel > 0)]
    acb = acb[~acb.allele.map(N.pg_label).str.endswith('P')]
    lone = acb.sort_values('count_panel', ascending=False).iloc[0]
    tnull = tm[null].sort_values('count', ascending=False).iloc[0]
    tma = pd.read_csv(C.TOMMO).query("gene == 'A'")
    harm = [('Genes compared', len(gl), 'harmonization.genes'),
            ('Controls with a genotype per gene, fewest', int(lt.n_samples.min()), 'harmonization.controls_min'),
            ('Controls with a genotype per gene, most', int(lt.n_samples.max()), 'harmonization.controls_max'),
            ('Our two-field alleles in controls (null alleles excluded)', int(cov.two_field.sum()),
             'harmonization.ours_two_field'),
            ('The same at P-group resolution', int(cov.ours_keys.sum()), 'harmonization.ours_groups'),
            ('ToMMo entries', len(tm), 'harmonization.tommo_entries'),
            ('ToMMo entries named by a P group', int(named_p.sum()), 'harmonization.tommo_named_p'),
            ('ToMMo entries for an allele in no P group', int((~named_p & ~null).sum()), 'harmonization.tommo_lone'),
            ('ToMMo null-allele entries', int(null.sum()), 'harmonization.tommo_null_entries'),
            ('Null-allele copies dropped, ours', int(lt.null_chr_ours.sum()), 'harmonization.null_copies_ours'),
            ('Null-allele copies dropped, ToMMo', int(lt.null_chr_panel.sum()), 'harmonization.null_copies_tommo'),
            ('Our alleles named after ToMMo\'s newest allele (P-group resolution)', int(cov.alleles_too_new.sum()),
             'harmonization.too_new_groups'),
            ('Newest IPD-IMGT/HLA release among ToMMo\'s alleles', f'{str(vint)[0]}.{str(vint)[1:3]}.{str(vint)[3]}',
             'tommo_vintage'),
            ('Example of merging: commonest allele named by another P group',
             f'{ex_merge} → {N.pg_label(ex_merge)} ({int(seen[ex_merge])} copies)', 'harmonization.example_merge'),
            ('Example of an allele in no P group, in both tables',
             f'{lone.allele} (ours {int(lone.count_ours)}, ToMMo {int(lone.count_panel)} copies)',
             'harmonization.example_lone'),
            ('Commonest ToMMo null entry', f'{tnull.allele} ({int(tnull["count"])} copies)', 'harmonization.example_null'),
            ('Rescaling example, ToMMo HLA-A commonest entry',
             f'{tma.sort_values("count").allele.iloc[-1]}: {int(tma["count"].max())} of {int(tma["count"].sum())} '
             f'copies, {int(tma.loc[~tma.allele.str.rstrip("P").str.endswith("N"), "count"].sum())} protein-coding',
             'harmonization.renorm')]
    harm = pd.DataFrame(harm, columns=['item', 'value', 'deck_key'])
    # the point off the diagonal (slide 19)
    out = freq_outlier(T, 'DRB1')
    # further checks, not on the slides: homozygosity / dropout, confidence against the panel, shift direction
    hz = T.hz.copy()
    cvp = T.cv.copy()
    sh = {'n': len(T.shift), 'kg_on_our_side': int((T.shift.kg_side == 'ours').sum()),
          'kg_on_tommo_side': int((T.shift.kg_side == 'ToMMo').sum())}
    unc = {'alleles': len(T.unc), 'copies': int(T.unc.count_ours.sum()), 'assessable_alleles': int(T.unc.assessable.sum()),
           'assessable_copies': int(T.unc[T.unc.assessable == 1].count_ours.sum())}
    return genes, harm, out, hz, cvp, sh, unc, vint


def freq_outlier(T, gene):
    """The typing deck's freq_outlier: the allele furthest from the diagonal (square-root axes) and why."""
    C, R = T.C, T.R
    al = T.ac[(T.ac.panel == 'ToMMo') & (T.ac.gene == gene)].copy()
    al['vis'] = (np.sqrt(al.freq_ours) - np.sqrt(al.freq_panel)).abs()
    o = al.sort_values('vis', ascending=False).iloc[0]
    hist = pd.read_csv(C.IMGT_HISTORY, comment='#', dtype=str).set_index('HLA_ID')
    stem = o.allele.split('*')[1]
    row = next((hist.loc[i] for i in hist.index if str(hist.loc[i].iloc[0]).startswith(f'{gene}*{stem}')
                and str(hist.loc[i].iloc[0]).rstrip('QNLSCA') == f'{gene}*{stem}'), None)
    out = [('Allele', o.allele), ('Our controls, frequency', float(o.freq_ours)), ('Our controls, copies', int(o.count_ours)),
           ('ToMMo, frequency', float(o.freq_panel)), ('ToMMo, copies', int(o.count_panel))]
    name = o.allele
    if row is not None:
        name = row.iloc[0]
        changed = next((rel for rel, v in row.items() if v != row.iloc[0]), None)
        out.append(('Name in IPD-IMGT/HLA 3.64.0', row.iloc[0]))
        if changed:
            rel = hist.columns[list(hist.columns).index(changed) - 1]
            out.append(('Renamed from', f'{row[changed]} (in release {rel[0]}.{rel[1:3]}.{rel[3]})'))
    pg = [ln.split(';') for ln in C.P_GROUPS.read_text().splitlines() if ln.startswith(f'{gene}*')]
    p_group = next((f[2] for f in pg if name.split('*')[1] in f[1].split('/') and f[2]), None)
    out.append(('P group', p_group or 'none'))
    freq = C.HLAHD_FREQ / f'{gene}_count.txt'
    out.append(('In HLA-HD\'s frequency table', 'yes' if (freq.exists() and f'{gene}*{stem}' in
                [w.split('\t')[0] for w in freq.read_text().split('\n')]) else 'no'))
    ref = pd.read_csv(R / '05.exports/residue_reference' / f'{gene}.tsv', sep='\t', index_col=0)
    if name in ref.index:
        a_ = ref.loc[name]
        bad = {'*', '.', '-', 'X'}
        best = None
        for other in ref.index:
            two = ':'.join(str(other).split(':')[:2])
            if two == name or two in (o.allele, name.rstrip('QNLSCA')):
                continue
            b_ = ref.loc[other]
            both = a_.notna() & b_.notna() & ~a_.isin(bad) & ~b_.isin(bad)
            n, d = int(both.sum()), int((a_[both] != b_[both]).sum())
            if n and (best is None or (d, -n) < (best[1], -best[2])):
                best = (two, d, n)
        if best and best[1] == 0:
            out.append(('Protein identical to', f'{best[0]} at all {best[2]} positions determined in both'))
    grp = al[al.allele.str.startswith(f'{gene}*{stem.split(":")[0]}:')]
    out.append((f'{gene}*{stem.split(":")[0]} group, our controls', float(grp.freq_ours.sum())))
    out.append((f'{gene}*{stem.split(":")[0]} group, ToMMo', float(grp.freq_panel.sum())))
    return pd.DataFrame(out, columns=['item', 'value'])


# ── T5 ───────────────────────────────────────────────────────────────────────
def t5(T):
    C, mk, summ = T.C, T.mk, T.summ
    genes = list(C.ASSOC_GENES)
    ad_cols, rd_cols = T.ad.columns, T.rd.columns
    enc = pd.DataFrame([{'gene': g, 'allele_columns': int((ad_cols.str.split('*').str[0] == g).sum()),
                         'residue_columns': int((rd_cols.str.split(':').str[0] == g).sum()),
                         'positions_polymorphic': int(T.sites.loc[T.sites.gene == g, 'polymorphic'].sum())}
                        for g in C.REPORTED_GENES])
    enc['assoc'] = enc.gene.isin(genes)
    # stepwise record per cohort x stage x gene (qc_record.tsv as step 06 wrote it)
    recs, per_gene = [], []
    for co in COHORTS:
        for key, step, name, rule in STAGES:
            r = T.qc[(T.qc.cohort == co) & (T.qc.step == step)]
            gl = C.REPORTED_GENES if key == 'typed' else genes
            for kind in ('allele', 'residue'):
                n = int(r.loc[(r['class'] == kind) & r.gene.isin(gl), 'columns'].sum())
                recs.append({'cohort': co, 'stage': key, 'stage_name': name, 'rule': rule, 'class': kind, 'columns': n})
        sc = summ['cohorts'][co]
        tot = {k: next(x['columns'] for x in recs if x['cohort'] == co and x['stage'] == 'maf' and x['class'] == k)
               for k in ('allele', 'residue')}
        if (tot['allele'], tot['residue']) != (sc['alleles_tested'], sc['residues_tested']):
            die(f'T5 {co}: the record and summary.json disagree')
    rec = pd.DataFrame(recs)
    # per gene: columns at each stage, copies covered (the deck's allele_coverage)
    ex = T.ad
    for co in COHORTS:
        cad = pd.read_csv(T.A / co / 'allele_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
        crd = pd.read_csv(T.A / co / 'residue_dosage.tsv', sep='\t', dtype={'sample_id': str}).set_index('sample_id')
        e = ex.loc[cad.index]
        if not (e[cad.columns].fillna(-1).astype(int) == cad.fillna(-1).astype(int)).all().all():
            die(f'T5 {co}: the cohort matrix is not a slice of the export')
        r = T.qc[T.qc.cohort == co]
        for g in C.REPORTED_GENES:
            row = {'cohort': co, 'gene': g, 'assoc': g in genes}
            for key, step, _, _ in STAGES:
                for kind in ('allele', 'residue'):
                    v = r.loc[(r.step == step) & (r['class'] == kind) & (r.gene == g), 'columns']
                    row[f'{kind}_{key}'] = int(v.sum()) if (key == 'typed' or g in genes) else np.nan
            if g in genes:
                gcols = [c for c in e.columns if c.split('*')[0] == g]
                typed = e[gcols].sum().sum()
                row['allele_copies_covered'] = e[[c for c in gcols if c in cad.columns]].sum().sum() / typed
            per_gene.append(row)
        # residue coverage at the kept positions
        er = T.rd.loc[crd.index]
        kept_pos = {c.rsplit(':', 1)[0] for c in crd.columns}
        at_kept = er[[c for c in er.columns if c.rsplit(':', 1)[0] in kept_pos]]
        per_gene.append({'cohort': co, 'gene': 'All', 'assoc': True,
                         **{f'{k}_{s}': int(rec[(rec.cohort == co) & (rec.stage == s) & (rec['class'] == k)]['columns'].iloc[0])
                            for s, _, _, _ in STAGES for k in ('allele', 'residue')},
                         'allele_copies_covered': sum(e[[c for c in e.columns if c.split('*')[0] == g and c in cad.columns]]
                                                      .sum().sum() for g in genes) /
                         sum(e[[c for c in e.columns if c.split('*')[0] == g]].sum().sum() for g in genes),
                         'residue_positions_kept': len(kept_pos),
                         'residue_copies_covered': float(crd.sum().sum() / at_kept.sum().sum())})
    pg = pd.DataFrame(per_gene)
    # floors (slides 22-23): call-rate and MAF distributions of the columns the floors act on
    edges = [0, 0.001, 0.005, 0.01, 0.05, 0.10, 0.25, 1.0]
    labels = ['<0.1%', '0.1-0.5%', '0.5-1%', '1-5%', '5-10%', '10-25%', '>=25%']
    fl = []
    for co in COHORTS:
        m = mk[mk.cohort == co]
        for kind in ('allele', 'residue'):
            k = m[(m['class'] == kind) & m.gene.isin(genes)]
            live = k[k.status != 'not_varying']
            cr = live[['call_rate_cases', 'call_rate_controls']].min(axis=1)
            ok = live[cr >= C.MIN_CALL_RATE]
            floors = {'mac20': ok.mac >= 20, 'maf05': ok.maf >= C.MIN_MAF, 'maf1': ok.maf >= 0.01}
            row = {'cohort': co, 'class': kind, 'columns': len(k), 'not_varying': int((k.status == 'not_varying').sum()),
                   'varying': len(live), 'pass_call_rate': len(ok), 'call_rate_min': float(cr.min()),
                   'call_rate_median': float(cr.median()),
                   'cr_lt90': int((cr < 0.90).sum()), 'cr_90_95': int(((cr >= 0.90) & (cr < 0.95)).sum()),
                   'cr_95_98': int(((cr >= 0.95) & (cr < 0.98)).sum()), 'cr_98_99': int(((cr >= 0.98) & (cr < 0.99)).sum()),
                   'cr_ge99': int((cr >= 0.99).sum()),
                   'lost_at_95': int((cr < 0.95).sum()), 'lost_at_98': int((cr < 0.98).sum()), 'lost_at_99': int((cr < 0.99).sum())}
            for b_, v in pd.cut(ok.maf, edges, right=False, labels=labels).value_counts().reindex(labels).items():
                row[f'maf_{b_}'] = int(v)
            for f, v in floors.items():
                row[f'kept_{f}'] = int(v.sum())
                row[f'floor_mac_{f}'] = int(ok.mac[v].min())
            fl.append(row)
    floors = pd.DataFrame(fl)
    # gene completeness (the ceiling of every marker's call rate) and the worst residue position
    ad = T.ad
    ph = T.pheno['full_mainland'].set_index('IID')
    comp = []
    rd = T.rd.loc[:, [c.split(':')[0] in genes for c in T.rd.columns]]
    pos = rd.columns.str.rsplit(':', n=1).str[0]
    worst = {}
    for pid in pd.unique(pos):
        cols = rd.columns[pos == pid]
        sm = rd[cols].sum(axis=1, min_count=len(cols))
        cr = float(((sm == 2).sum()) / len(rd))
        g = pid.split(':')[0]
        if g not in worst or cr < worst[g][0]:
            worst[g] = (cr, pid.split(':')[1])
    d_full = ad.loc[ph.index]
    for g in genes:
        cols = [c for c in ad.columns if c.split('*')[0] == g]
        has_all = ad[cols].sum(axis=1, min_count=1).notna()
        has = d_full[cols].sum(axis=1, min_count=1).notna()
        row = {'gene': g, 'all_genomes': float(has_all.mean()), 'cases_full': float(has[ph.PHENO1 == 1].mean()),
               'controls_full': float(has[ph.PHENO1 == 0].mean()), 'worst_position': worst[g][1] if g in worst else '',
               'worst_position_call_rate': worst[g][0] if g in worst else np.nan}
        for co in COHORTS:
            ids = [i for i in T.pheno[co].IID if i in ad.index]
            row[co] = float(ad.loc[ids, cols].sum(axis=1, min_count=1).notna().mean())
        comp.append(row)
    comp = pd.DataFrame(comp)
    # residue positions losing call rate to the alignment (full mainland)
    adf, rdf = ad.loc[ph.index], T.rd.loc[ph.index]
    loss = []
    for g in genes:
        cg = adf[[c for c in adf.columns if c.split('*')[0] == g]].notna().all(axis=1).mean()
        cols = [c for c in rdf.columns if c.split(':')[0] == g]
        for _, cc in pd.Series(cols).groupby(pd.Series(cols).map(lambda c: c.rsplit(':', 1)[0])):
            loss.append(cg - rdf[list(cc)].notna().all(axis=1).mean())
    loss = pd.Series(loss)
    extra = {'residue_positions': len(loss), 'positions_losing_over_1pct': int((loss > 0.01).sum()),
             'positions_losing_over_5pct': int((loss > 0.05).sum()),
             'polymorphic_positions_assoc_genes': int(T.sites[T.sites.gene.isin(genes)].polymorphic.sum()),
             'allele_columns_exported': int(len(ad_cols)), 'residue_columns_exported': int(len(rd_cols)),
             'expected_case_copies': {co: {f: round(f * 2 * summ['cohorts'][co]['cases'], 1) for f in (0.005, 0.01)}
                                      for co in COHORTS},
             'cohorts': {co: {k: summ['cohorts'][co][k] for k in ('samples', 'cases', 'controls', 'alleles_tested',
                                                                  'residues_tested')} for co in COHORTS},
             'thresholds': summ.get('thresholds')}
    return enc, rec, pg, floors, comp, extra


# ── T6 ───────────────────────────────────────────────────────────────────────
def t6(T, compared):
    C, N = T.C, T.N
    pm = T.pgmap.set_index('allele').p_group
    splits, rows = [], []
    X = T.ad.loc[T.controls_full()]
    acx = T.ac[T.ac.panel == 'ToMMo'].set_index(['gene', 'allele'])
    for co in COHORTS:
        m = T.mk[(T.mk.cohort == co) & (T.mk['class'] == 'allele') & (T.mk.status == 'tested')].copy()
        m['pg'] = m.marker.map(pm)
        if m.pg.isna().any():
            die(f'T6 {co}: tested allele without a P group')
        size = m.groupby('pg').marker.transform('size')
        sel = m[size > 1]
        cmp_sel = sel[sel.gene.isin(compared)]
        splits.append({'cohort': co, 'tested': len(m), 'groups': int(sel.pg.nunique()), 'columns': len(sel),
                       'groups_compared': int(cmp_sel.pg.nunique()), 'columns_compared': len(cmp_sel)})
        if co != 'full_mainland':
            continue
        sel = sel.assign(freq=sel.copies / (2 * sel.called), rank=sel.gene.map(list(C.ASSOC_GENES).index))
        for pg_, d in sel.sort_values(['rank', 'pg']).groupby('pg', sort=False):
            g = d.gene.iloc[0]
            rep = pg_[:-1] if pg_.endswith('P') else pg_
            d = d.assign(first=d.marker != rep).sort_values(['first', 'freq'], ascending=[True, False])
            cols = [c for c in X.columns if c.split('*')[0] == g and not c.endswith('N')]
            s = X[cols].sum()
            s, tot = s[s > 0], float(s.sum())
            key = N.pg_of(d.marker.iloc[0])
            members = s[[N.pg_of(a) == key for a in s.index]]
            tested = list(d.marker)
            others = members.drop([a for a in tested if a in members.index])
            row = {'gene': g, 'p_group': pg_, 'compared': g in compared,
                   'tested': '; '.join(f'{a}' for a in tested),
                   'tested_freq_controls': '; '.join(f'{s.get(a, 0) / tot:.6f}' for a in tested),
                   'tested_freq_all': '; '.join(f'{f:.6f}' for f in d.freq),
                   'others_n': int(len(others)), 'others_freq': float(others.sum() / tot),
                   'group_freq_controls': float(members.sum() / tot),
                   'tommo_freq': float(acx.loc[(g, key), 'freq_panel']) if g in compared else np.nan}
            rows.append(row)
    # the DQA1 worked example (slide 29)
    gene = 'DQA1'
    a = T.ac[(T.ac.panel == 'ToMMo') & (T.ac.gene == gene)]
    s = X[[c for c in X.columns if c.split('*')[0] == gene and not c.endswith('N')]].sum()
    s = s[s > 0]
    members = {}
    for al_, n in s.items():
        members.setdefault(N.pg_of(al_), []).append((al_, int(n)))
    ex = []
    for _, r in a.iterrows():
        where = 'both' if r.count_ours > 0 and r.count_panel > 0 else ('ours' if r.count_ours > 0 else 'tommo')
        mm = sorted(members.get(r.allele, []), key=lambda x: -x[1])
        if sum(n for _, n in mm) != int(r.count_ours):
            die(f'T6 {r.allele}: members do not add up')
        ex.append({'p_group': N.pg_label(r.allele), 'where': where, 'freq_ours': float(r.freq_ours),
                   'freq_tommo': float(r.freq_panel), 'copies_ours': int(r.count_ours), 'copies_tommo': int(r.count_panel),
                   'members': ', '.join(f'{x} ({n})' for x, n in mm)})
    ex = pd.DataFrame(ex)
    ex['order'] = ex['where'].map({'both': 0, 'ours': 1, 'tommo': 2})
    ex = ex.sort_values(['order', 'freq_tommo', 'freq_ours'], ascending=[True, False, False]).drop(columns='order')
    return pd.DataFrame(splits), pd.DataFrame(rows), ex


# ── T7 / T8 ──────────────────────────────────────────────────────────────────
STATUS = {'tested': 'Tested', 'maf': 'Below MAF 0.5%', 'call_rate': 'Call rate below 95%',
          'excluded_gene': 'Gene not tested', 'not_varying': 'Not varying'}


def t7_t8(T):
    C, N = T.C, T.N
    mk = T.mk.copy()
    out = {}
    for kind in ('allele', 'residue'):
        m = mk[mk['class'] == kind]
        base = m[m.cohort == 'full_mainland'][['marker', 'gene', 'position', 'residue']].copy()
        for co in COHORTS:
            d = m[m.cohort == co].set_index('marker')
            ph = T.pheno[co].set_index('IID')
            mat = (T.ad if kind == 'allele' else T.rd).loc[ph.index]
            case = mat[ph.PHENO1.values == 1]
            ctrl = mat[ph.PHENO1.values == 0]
            c = COHORTS.index(co)
            tag = ['n', 'i', 'f'][c]
            base[f'{tag}_status'] = base.marker.map(d.status).map(STATUS)
            base[f'{tag}_cr_cases'] = base.marker.map(d.call_rate_cases)
            base[f'{tag}_cr_controls'] = base.marker.map(d.call_rate_controls)
            base[f'{tag}_copies_cases'] = base.marker.map(case.sum())
            base[f'{tag}_copies_controls'] = base.marker.map(ctrl.sum())
            base[f'{tag}_maf'] = base.marker.map(d.maf)
            if (base[f'{tag}_copies_cases'] + base[f'{tag}_copies_controls'] != base.marker.map(d.copies)).any():
                die(f'T7/T8 {co}: copies by group do not add up to the record')
        if kind == 'allele':
            pg = T.pgmap.set_index('allele')
            mt = T.match.set_index('allele')
            base['p_group'] = base.marker.map(lambda a: N.pg_label(a))
            base['first_release'] = base.marker.map(lambda a: N.release_of(a))
            base['match'] = base.marker.map(mt.how).map(lambda h: 'exact' if h == 'exact' else
                                                        ('deeper allele' if str(h).startswith('prefix') else 'none'))
            base['matched_to'] = base.marker.map(mt.matched_to)
            ton = T.ac[T.ac.panel == 'ToMMo'].set_index(['gene', 'allele'])
            kg = T.ac[T.ac.panel == '1KG JPT'].set_index(['gene', 'allele'])
            key = base.marker.map(N.pg_of)
            base['tommo_freq'] = [float(ton.loc[(g, k), 'freq_panel']) if (g, k) in ton.index else np.nan
                                  for g, k in zip(base.gene, key)]
            base['tommo_compared'] = base.gene.isin(set(T.loc[(T.loc.panel == 'ToMMo') &
                                                                (T.loc.frequency_comparable == 1)].gene))
            base['kg_freq'] = [float(kg.loc[(g, k), 'freq_panel']) if (g, k) in kg.index else np.nan
                               for g, k in zip(base.gene, key)]
            base['kg_compared'] = base.gene.isin(set(T.loc[T.loc.panel == '1KG JPT'].gene))
            if len(base) != T.ad.shape[1]:
                die('T7: not every exported allele column has a row')
        else:
            if len(base) != T.rd.shape[1]:
                die('T8: not every exported residue column has a row')
        out[kind] = base
    return out['allele'], out['residue']


def main():
    a = parse_args()
    o = Path(a.out)
    T = Typing(a)
    run, plat, reasons = t1(T)
    write_tsv(run, o / 't1_run.tsv')
    write_tsv(plat, o / 't1_platforms.tsv')
    drb, two = t2(T)
    write_tsv(drb, o / 't2_drb.tsv')
    write_tsv(two, o / 't2_two_allele.tsv')
    genes, tot = t3(T)
    write_tsv(genes, o / 't3_genes.tsv')
    t4g, harm, outl, hz, cvp, sh, unc, vint = t4(T)
    write_tsv(t4g, o / 't4_panels.tsv')
    write_tsv(harm, o / 't4_harmonization.tsv')
    write_tsv(outl, o / 't4_outlier.tsv')
    write_tsv(hz, o / 't4_homozygosity.tsv')
    write_tsv(cvp, o / 't4_confidence_vs_panel.tsv')
    enc, rec, pg, floors, comp, extra = t5(T)
    write_tsv(enc, o / 't5_encoding.tsv')
    write_tsv(rec, o / 't5_record.tsv')
    write_tsv(pg, o / 't5_genes.tsv')
    write_tsv(floors, o / 't5_floors.tsv')
    write_tsv(comp, o / 't5_completeness.tsv')
    compared = list(t4g.loc[(t4g.panel == 'ToMMo 61KJPN') & (t4g.frequency_comparable == 1), 'gene'])
    sp, rows, ex = t6(T, compared)
    write_tsv(sp, o / 't6_split_counts.tsv')
    write_tsv(rows, o / 't6_split_groups.tsv')
    write_tsv(ex, o / 't6_dqa1.tsv')
    t7, t8 = t7_t8(T)
    write_tsv(t7, o / 't7_alleles.tsv')
    write_tsv(t8, o / 't8_residues.tsv')
    write_json({'totals': tot, 'outlier_reasons': reasons, 'shift': sh, 'absent_from_tommo': unc,
                'tommo_vintage': vint, 'markers': extra, 'compared_genes': compared,
                'reported_genes': list(T.C.REPORTED_GENES), 'assoc_genes': list(T.C.ASSOC_GENES)},
               o / 'typing_facts.json')
    log(f'typing tables written: T7 {len(t7)} alleles, T8 {len(t8)} residues')


if __name__ == '__main__':
    main()
