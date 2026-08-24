#!/usr/bin/env python3
# ==============================================================================
# figure_catalogue.py — ONE document per figure family per cohort.
#
# Why:      The per-locus figures are produced by a fan-out, so a scan with 100
#           peaks produced 100 regional figures, 100 fine-mapping figures and 100
#           conditional figures — and, when every figure wrote its own sidecar,
#           300 documents whose prose was identical and whose only difference was
#           a table of six to thirteen numbers. That is not documentation; it is
#           the same page copied 300 times, and it buries the handful of
#           documents a reader actually needs.
#
#           So the prose is written ONCE per family, and the per-locus numbers
#           become rows of one table inside it. A reader opens one file per
#           family and can compare loci down a column, which the per-locus files
#           made impossible.
#
# Input:    the `*.stats.json` each locus figure drops beside its PNG — the same
#           `numbers` list those figures used to render into their own sidecar.
# Output:   <figures>/<NN.family>/README.md
#
# The prose lives here rather than in the plotting scripts because it describes
# the family, not any one rendering. A plotting script draws; this explains.
# ==============================================================================
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

import figure_doc
import plot_style as S

# Peak ids are `gw###_<chrom>_<pos>` / `sg###_<chrom>_<pos>`; the prefix is the
# tier and is load-bearing across both pipelines (see call_peaks.py).
TIERS = {'sig': 'significant', 'sug': 'suggestive',   # genes, tiered by BH
         'gw': 'genome_wide', 'sg': 'suggestive'}     # variant peaks, fixed P


# Families whose members carry no tier of their own. Their ids are scan names,
# not gene or peak ids, so tier_of() has nothing to read.
UNTIERED = {'scan'}


def tier_of(peak_id):
    s = str(peak_id)
    for n in (3, 2):                 # longest first: 'sig' must not read as 'sg'
        t = TIERS.get(s[:n])
        if t is not None:
            return t
    raise SystemExit(f"cannot determine the tier of '{peak_id}' — expected the "
                     'gw/sg prefix call_peaks.py writes or the sig/sug prefix '
                     'gene_scan.py writes')


# ── family prose ─────────────────────────────────────────────────────────────
# Each entry returns the keyword arguments figure_doc.write_doc takes. `args`
# carries the run's own configuration so nothing here hardcodes a threshold, a
# covariate set or a panel's data source.

def _regional(args):
    titles = dict(S.LD_SOURCE_TITLE,
                  **dict(s.split('=', 1) for s in args.panel_label if '=' in s))
    return dict(
        title=f'Regional association — {args.cohort}',
        question=('Is the LD structure that produces each peak a property of the East Asian '
                  'population, or an artefact of this particular sample?'),
        panels=[
            ('a', 'In-sample cohort LD',
             'The association statistics of the peak window, each variant coloured by r2 with the '
             'lead computed on the *same* samples the statistics came from. This is the LD matrix the '
             'SuSiE fine-mapping uses, so this panel shows exactly what fine-mapping saw.'),
            ('b', titles['tommo'],
             f'Identical statistics, recoloured by r2 from the {titles["tommo"]} co-occurrence '
             'tables. That resource publishes only pairs with r2 >= 0.2, so a variant present in the '
             'panel but returning no pair with the lead is *bounded* below 0.2 — it belongs in the '
             'lowest colour bin, not in grey. Grey is reserved for variants absent from the panel '
             'entirely.'),
            ('c', titles['1000g_eas'],
             f'Identical statistics again, coloured by r2 in {titles["1000g_eas"]} — an '
             'out-of-sample population reference with no relationship to this cohort.'),
            ('d', 'Gene models',
             'Filled boxes are exons, the connecting line spans the introns and chevrons give the '
             'transcribed strand. One representative transcript per gene (the longest protein-coding '
             'one), because drawing every transcript of a multi-transcript gene would need many rows '
             'to say the same thing. Only genes with an official symbol appear; clone-accession '
             'models name a sequencing clone rather than a gene and are suppressed.'),
        ],
        reading=[
            'Compare the three colour patterns. If they agree, the LD block is a population property '
            'and the credible set can be trusted to reflect real correlation structure.',
            'If the in-sample panel alone shows tight LD, the correlation is being driven by sampling '
            'noise at this N and the credible set is correspondingly fragile.',
            'Read the lead against the gene track: whether it sits in an exon, an intron or between '
            'genes constrains which mechanisms are plausible.',
            'Both tier lines are drawn in every panel, so a suggestive locus can be read against the '
            'genome-wide line it did not reach.',
        ],
        limits=[
            'It does not identify a causal variant. LD colour is correlation with the lead, not '
            'evidence of function — that is what the fine-mapping figures address.',
            'The two external panels differ from each other and from the study samples in ancestry '
            'breadth and in genomic coverage. Disagreement between (b) and (c) may reflect that '
            'difference rather than an error in either.',
        ],
        defs=['ld_r2', 'gw_sig'],
        interpretation=('Agreement between the three panels means the LD structure is a property of '
                        'the population rather than of this sample; divergence localised to the '
                        'in-sample panel would indicate that the credible set is being shaped by '
                        'sampling noise at this N. Each panel prints its own coverage in the table '
                        'below, so a sparse source cannot be mistaken for a low-LD region.'))


def _finemap(args):
    return dict(
        title=f'Fine-mapping — {args.cohort}',
        question='Does the posterior concentrate on a few variants, or does LD spread it out?',
        panels=[
            ('a', 'Association',
             'The peak window, coloured by in-sample r2 with the lead (purple diamond). Shown so the '
             'PIP panel can be read against the evidence that produced it. Both tier lines are drawn.'),
            ('b', 'Posterior inclusion probability',
             'Per-variant PIP from susie_rss on the summary statistics and the in-sample signed-r LD '
             'matrix. Rings mark credible-set membership; the bracket spans a set\'s physical extent. '
             'A credible set is the smallest group of variants carrying 95% of the posterior mass for '
             'one signal.'),
            ('c', 'Resolution',
             'Cumulative posterior mass against the PIP-ranked members of each set. A curve reaching '
             '0.95 in a few steps means the signal is resolved to those variants; a slow curve means '
             'LD has spread the mass and the set cannot be narrowed at this sample size.'),
        ],
        reading=[
            'A credible set of one or two variants is a resolved signal. A set of twenty means the '
            'data cannot distinguish among them, not that twenty variants are causal.',
            'Check that the set contains the lead. If it does not, the lead is a tag and the '
            'posterior prefers a neighbour.',
            'PIP is conditional on the LD matrix being the one the statistics came from, which is why '
            'in-sample LD is used here rather than a reference panel.',
            'Read `credible sets` = 0 in the table as "SuSiE found no set meeting the purity and '
            'coverage requirements", not as "no signal".',
        ],
        limits=[
            'SuSiE assumes the causal variant is present in the data. A causal variant not genotyped '
            'or filtered out cannot appear, and its posterior mass will be distributed over its tags.',
            'It cannot rank the biological plausibility of set members — only their statistical '
            'compatibility with the observed association pattern.',
        ],
        defs=['pip', 'ld_r2'],
        interpretation=('Fine-mapping is conditional on the LD matrix being the one the statistics '
                        'came from, which is why in-sample LD is used here. At this effective sample '
                        'size the posterior is driven by a small number of strongly associated '
                        'variants, so a wide credible set should be read as insufficient resolution '
                        'rather than as evidence against a single causal variant.'),
        model=S.FORMULAS['susie'])


def _conditional(args):
    return dict(
        title=f'Conditional analysis — {args.cohort}',
        question='Does each peak carry one association signal, or more than one?',
        panels=[
            ('a', 'The peak, by conditioning round',
             'Round 0 is the unconditioned fit restricted to the peak window. Each later round adds '
             'the previous round\'s top variant to the covariate set and re-fits every variant in the '
             'window. Same samples, same model and the same covariate set as the genome-wide scan, so '
             'the rounds are comparable to it and to each other.'),
            ('b', 'The stepwise decision',
             'The statistic the procedure actually acted on: the smallest P remaining at each round, '
             'against the threshold of the peak\'s own tier. The loop stops at the first round in '
             'which nothing in the window clears it.'),
        ],
        reading=[
            'If the whole peak collapses after conditioning on the lead, every significant variant '
            'there was tagging one underlying signal — the expected result inside a single LD block.',
            'A residual peak that still clears the threshold is a second, independent signal, and the '
            'procedure will have added it to the conditioning set and continued.',
            'Panel (b) is the audit trail: it shows the number the stopping rule compared, round by '
            'round.',
            f'Both tier lines are drawn; the one the peak was judged against — '
            f'{S.sci(args.alpha)} for a genome-wide peak, {S.sci(args.alpha_suggestive)} for a '
            'suggestive one — is drawn heavier and marked "decision".',
        ],
        limits=[
            'It cannot separate two causal variants in near-perfect LD; conditioning on one removes '
            'both. That limit is a property of the sample, not of the method.',
            'Absence of a secondary signal at this N is weak evidence — the conditional test faces the '
            'same detection floor as the primary scan, applied to a residual effect.',
            'A suggestive peak is judged against the suggestive threshold, so its signal count is not '
            'comparable to a genome-wide peak\'s without accounting for the different bar.',
        ],
        defs=['conditional', 'gw_sig', 'model'],
        interpretation=('Collapse of the whole locus after conditioning on the lead means every '
                        'significant variant there was tagging one underlying signal — the expected '
                        'outcome inside a single LD block. A residual peak would mark a second, '
                        'independent signal. Conditioning is performed on the same samples, model and '
                        'covariates as the genome-wide scan, so the rounds are comparable to it and '
                        'to each other.'),
        model=S.glm_formula(args.pc_label, args.n_pcs))


def _gene(args):
    return dict(
        title=f'Gene-based hits — {args.cohort}',
        question=('For each gene the burden test flagged, how many carriers actually '
                  'produced the signal, and is the difference between cases and controls '
                  'or between the cohort and the population?'),
        panels=[
            ('a', 'Case against control, per variant',
             'Alternate-allele frequency of every variant the gene-based test collapsed, in '
             'cases and in controls. A gene whose signal rests on ONE variant is a different '
             'claim from one where several variants agree, and the gene-level P cannot '
             'distinguish them.'),
            ('b', 'Control against the population reference',
             'The same variants, control frequency against ToMMo 60KJPN. This is the artefact '
             'check that matters under this design: cases and controls were sequenced on '
             'different platforms at different depth, so a control frequency that departs from '
             'the population reference is a genotyping difference rather than a disease '
             'association. Agreement here is what makes the case excess in (a) worth reading.'),
            ('c', 'The collapsed burden',
             'Cumulative minor-allele frequency over the gene in each group, labelled with the '
             'allele counts the test actually saw.'),
        ],
        reading=[
            'Read the allele counts in (c) first. A gene-level P of 1e-7 computed from a dozen '
            'carriers is a statement about a dozen people, and the odds ratio that comes with '
            'it has an interval spanning an order of magnitude.',
            'Check (b) before believing (a). Under a fully confounded design the control side '
            'is where a technical artefact shows, and ToMMo is the only external reference '
            'available.',
            'In (a), see whether the excess is spread across variants or carried by one. '
            'Burden tests weight every qualifying variant equally, so a single variant can '
            'carry a whole gene.',
        ],
        limits=[
            'It cannot establish that a variant is causal, or that the gene is. A burden test '
            'reports that qualifying variants are unequally distributed between groups.',
            'ToMMo agreement rules out one artefact class, not all of them. Residual platform '
            'chemistry is not fixable from within this data set.',
            'The odds ratio is fitted by logistic regression on the same covariates as the '
            'scan, but at these allele counts it is a description of the carriers observed, '
            'not a population estimate.',
        ],
        defs=['or', 'eaf'],
        interpretation=('A gene-based hit under this design is a candidate, not a confirmed '
                        'association. tuning.rv fixed the minor-allele-count threshold at the '
                        'point where the adjusted apparent group effect reaches a statistical '
                        'null, which removes the genome-wide burden artefact; it does not make '
                        'any individual gene real. What these figures add to the P-value is the '
                        'number of people behind it and an external frequency check.'))


def _scan(args):
    return dict(
        title=f'Gene-based scans — {args.cohort}',
        question=('Is each scan calibrated, and which genes does it single out?'),
        panels=[
            ('a', 'Gene-based association across the genome',
             'One point per gene, at its midpoint, banded by chromosome. The two lines are '
             'the P-values that Benjamini-Hochberg implies for FDR 0.05 and 0.10 IN THIS '
             'SCAN. They are not shared constants: BH adapts to the number of genes tested, '
             'and that ranges from ~160 under HIGH to ~13,000 under LOW+MODERATE+HIGH, so '
             'the implied P differs by roughly an order of magnitude between strata. Each '
             'figure therefore prints its own.'),
            ('b', 'Quantile-quantile plot',
             'The same gene P-values against a uniform null, with a Beta(i, n-i+1) 95% band '
             'and the gene-level genomic-control inflation factor. Read over GENES, not '
             'variants — it measures whether the gene-level null is calibrated.'),
        ],
        reading=[
            'Read lambda first. A scan whose gene-level lambda is far from 1 has an '
            'inflated null, and its tier counts are correspondingly optimistic.',
            'Compare strata down the table below rather than across the figures: the same '
            'tier label means a different P in each stratum, and the table states both.',
            'A gene named in red cleared FDR 0.05; grey cleared only 0.10.',
        ],
        limits=[
            'It cannot separate real polygenic signal from residual confounding. Under this '
            'design cases and controls differ in sequencing platform and depth, and the '
            'minAC inherited from tuning.rv removes the genome-wide burden artefact, not '
            'every per-gene one.',
            'A tier is not a finding. There is no independent cohort for this phenotype.',
        ],
        defs=['lambda_gc'],
        interpretation=('The scans are the shape of the data, not a list of results. What a '
                        'reader should take from them is the calibration in (b) and the '
                        'handful of genes that stand clear of the BH line in (a); the '
                        'per-gene figures are where the carriers behind each of those '
                        'become visible.'))


FAMILIES = {'regional': _regional, 'finemap': _finemap,
            'conditional': _conditional, 'gene': _gene, 'scan': _scan}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--family', required=True, choices=sorted(FAMILIES))
    p.add_argument('--cohort', required=True)
    p.add_argument('--stats-dir', required=True,
                   help='Directory holding the per-locus *.stats.json files.')
    p.add_argument('--out', required=True)
    p.add_argument('--figure-dir', required=True,
                   help='Path the figures are published under, quoted in the header.')
    p.add_argument('--alpha', type=float, default=S.GW_ALPHA)
    p.add_argument('--alpha-suggestive', type=float, default=S.SUGGEST_ALPHA)
    p.add_argument('--pc-label', default='')
    p.add_argument('--n-pcs', type=int, default=0)
    p.add_argument('--panel-label', action='append', default=[])
    p.add_argument('--methods-ref', default='../../../../docs/METHODS.md')
    return p.parse_args()


def main():
    args = parse_args()
    files = sorted(Path(args.stats_dir).glob(f'{args.family}.*.stats.json'))
    rows = []
    for f in files:
        d = json.loads(f.read_text())
        row = {'peak': d['peak']}
        if args.family not in UNTIERED:
            row['tier'] = tier_of(d['peak'])
        row.update({k: v for k, v in d['values']})
        rows.append(row)

    # Sorted so the genome-wide loci lead and, within a tier, the peak ids run in
    # the order call_peaks assigned them — the same order every other table uses.
    df = pd.DataFrame(rows)
    if len(df) and 'tier' not in df.columns:
        df = df.sort_values('peak')
        n_gw = n_sg = 0
    elif len(df):
        df = df.sort_values(['tier', 'peak'], key=lambda c: c.map(
            {'genome_wide': 0, 'suggestive': 1}).fillna(2) if c.name == 'tier' else c)
        # `cohort` repeats on every row and is already in the title.
        df = df.drop(columns=[c for c in ('cohort',) if c in df.columns])
        n_gw = int((df['tier'] == 'genome_wide').sum())
        n_sg = int((df['tier'] == 'suggestive').sum())
    else:
        n_gw = n_sg = 0

    kw = FAMILIES[args.family](args)
    tier_path = '' if args.family in UNTIERED else '<tier>/'
    subject = (f'**Figures:** `{args.figure_dir}/{tier_path}{args.family}.<peak>.png` — '
               f'{len(df)} locus figure(s): {n_gw} genome-wide, {n_sg} suggestive.\n\n'
               'One document covers the whole family: the panels, the reading and the limits are '
               'properties of the figure, identical for every locus. Only the numbers differ, and '
               'they are tabulated below, one row per locus.')
    figure_doc.write_doc(
        None, out_path=args.out, subject=subject,
        numbers=[('cohort', args.cohort), ('loci', len(df)),
                 ('genome-wide loci', n_gw), ('suggestive loci', n_sg),
                 ('genome-wide threshold', S.sci(args.alpha)),
                 ('suggestive threshold', S.sci(args.alpha_suggestive))],
        tables=[('Per-locus values', df)] if len(df) else None,
        methods_ref=args.methods_ref, **kw)
    print(f'[figure_catalogue] {args.cohort} {args.family}: '
          f'{len(df)} loci ({n_gw} genome-wide, {n_sg} suggestive) -> {args.out}')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error in figure_catalogue: {e}', file=sys.stderr)
        sys.exit(1)
